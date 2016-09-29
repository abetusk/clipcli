#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <unistd.h>
#include <string.h>

#include <errno.h>

#include <vector>

#include "clipper.hpp"
#include "poly2tri.h"

using namespace std;
using namespace ClipperLib;

//bool g_verbose_flag=false;
int g_verbose_level = 0;
char gVersion[]="0.0.1";

char *g_function=NULL;

long long int gMulFactor=1;
double gOffsetMiterLimit = 3.0;
double gOffsetRadius = 0.0;
double gEps = 0.000001;
int gForceOrientation = 0;

int gReadFloatFlag = 0;

int gSortOrder = 0;
int gPolyTreeFlag = 0;

bool gPrintReverse = false;
bool gPrintBoundingBox = false;
int gBoundingBoxMargin = 0;

void show_help(void) {
  printf("usage:\n\n");
  printf("clipcli [-t clip-type] [-S subj-fill] [-C clip-fill] [-v] [-V] [-s subj] [-s subj] ... [-c clip] [-c clip] ...\n\n");
  printf("  [-t clip-type]      Clipping operation.  One of (Intersection, Union, Difference, Xor).  Union default.\n");
  printf("  [-S subj-fill]      Subject fill type.  One of (EvenOdd, NonZero, Positive, Negative).  NonZero default.\n");
  printf("  [-C clip-fill]      Clip fill type.  One of (EvenOdd, NonZero, Positive, Negative).  NonZero default.\n");
  printf("  [-s subj]           File containg subject polygon.\n");
  printf("  [-c clip]           File containg clip polygon.\n");
  printf("  [-x mul]            Apply multiplication factor to input polygons.\n");
  printf("  [-R radius]         Polygon offset radius.\n");
  printf("  [-M miter_limit]    Miter limit (default %f).\n", gOffsetMiterLimit);
  printf("  [-E epsilon]        Epsilon (default %f).\n", gEps);
  printf("  [-O orientation]    force orientation of all polygons (>0 cw, <0 ccw, 0 (default) input orientation).\n");
  printf("  [-F]                read floating point input instead of long long integer\n");
  printf("  [-v]                Increase verbose level (can be specified multiple times).\n");
  printf("  [-P sortorder]      Specify what order to output polygons.  One of (area-inc, area-dec, cw-inc, cw-dec, pre, post, undef, infix).  Default undef/infix.\n");
  printf("  [-B]                Print bounding box (cw)\n");
  printf("  [-b margin]         Add margin to bounding box\n");
  printf("  [-r]                print reverse order polygon contours\n");
  printf("  [-T]                Output in tree order\n");
  printf("  [-V]                Show version.\n");
  printf("\n");

}

void processSortOrder(char *so) {
  gSortOrder = 0;

  if (strcmp(so, "area-inc")==0) {
    gSortOrder = 1;
  }
  else if (strcmp(so, "area-dec")==0) {
    gSortOrder = 2;
  }
  else if ((strcmp(so, "cw-inc")==0) || (strcmp(so,"ccw-dec")==0)) {
    gSortOrder = 3;
  }
  else if ((strcmp(so, "cw-dec")==0) || (strcmp(so,"ccw-inc")==0))  {
    gSortOrder = 4;
  } else if (strcmp(so, "pre")==0) {
    gSortOrder = -1;
  } else if (strcmp(so, "post") == 0) {
    gSortOrder = 1;
  }
}

typedef struct sort_order_type {
  int orig_idx;
  double val;
  unsigned long long int X, Y;
} sort_order_t;

int sort_order_cmp(const void *x, const void *y) {
  int inc_dec = 1;
  sort_order_t *a, *b;
  a = (sort_order_t *)x;
  b = (sort_order_t *)y;

  if ((gSortOrder==1) || (gSortOrder==3)) { inc_dec = 1; }
  else { inc_dec = -1; }

  if ((a->val) == (b->val)) { return 0; }
  if ((a->val) < (b->val)) { return inc_dec; }
  return -inc_dec;
}

int sort_order_xy_cmp(const void *x, const void *y) {
  int inc_dec = 1;
  sort_order_t *a, *b;
  a = (sort_order_t *)x;
  b = (sort_order_t *)y;

  if ((gSortOrder==1) || (gSortOrder==3)) { inc_dec = 1; }
  else { inc_dec = -1; }

  if (((a->X) == (b->X)) && ((a->Y) == (b->Y))) { return 0; }

  if (a->X < b->X) { return inc_dec; }
  else if (a->X > b->X) { return -inc_dec; }

  if (a->Y < b->Y) { return inc_dec; }
  else if (a->Y > b->Y) { return -inc_dec; }

  return 0;
}


void show_version(void) {
  printf("%s\n", gVersion);
}

void load_poly( FILE *fp, Paths &p ) {
  int i, j, k;
  int line_no=0;
  long long int X, Y;
  char buf[1024];
  bool first=true;
  bool is_outer_boundary = true;
  double X_d, Y_d;

  Path cur_clip_path, clip_apths_res;
  
  while (fgets(buf, 1024, fp)) {
    line_no++;

    if (buf[0] == '#') continue;
    if (buf[0] == '\n') {
      if (cur_clip_path.size() == 0) { continue; }
      p.push_back( cur_clip_path );
      cur_clip_path.clear();
      continue;
    }


    if (gReadFloatFlag) {
      k = sscanf(buf, "%lf%*[ \t]%lf", &X_d, &Y_d);
      if (k!=2) { perror(""); fprintf(stderr, "invalid read on line %i, exiting\n", line_no); exit(2); }
      X = (gMulFactor)*X_d;
      Y = (gMulFactor)*Y_d;
      cur_clip_path.push_back( ClipperLib::IntPoint(X, Y) );
    } else {
      k = sscanf(buf, "%lli%*[ \t]%lli", &X, &Y);
      if (k!=2) { perror(""); fprintf(stderr, "invalid read on line %i, exiting\n", line_no); exit(2); }
      cur_clip_path.push_back( ClipperLib::IntPoint(X*gMulFactor, Y*gMulFactor) );
    }


  }

  if (cur_clip_path.size() > 0) {

    if (gForceOrientation!=0) {
      if (gForceOrientation>0) {
        if (!ClipperLib::Orientation(cur_clip_path)) {
          ClipperLib::ReversePath(cur_clip_path);
        }
      } else {
        if (ClipperLib::Orientation(cur_clip_path)) {
          ClipperLib::ReversePath(cur_clip_path);
        }
      }
    }

    p.push_back(cur_clip_path);

  }

}


void load_polys( vector<char *> &fns, Paths &polys) {
  FILE *fp;
  int i, j, k;

  for (i=0; i<fns.size(); i++) {
    fp = fopen(fns[i], "r");
    load_poly( fp, polys );
    fclose(fp);
  }
}

ClipType find_clip_type( char *f ) {
  int i, j, k;
  char x[][16] = { "Intersection", "Union", "Difference", "Xor",
                   "intersection", "union", "difference", "xor",
                   "int", "un", "diff", "xor",
                   "i", "u", "d", "x" };
  ClipType lookup[16] = { ctIntersection, ctUnion, ctDifference, ctXor,
                          ctIntersection, ctUnion, ctDifference, ctXor,
                          ctIntersection, ctUnion, ctDifference, ctXor,
                          ctIntersection, ctUnion, ctDifference, ctXor };
  for (i=0; i<16; i++) {
    if (strncmp( x[i], f, 16 ) == 0) {
      return lookup[i];
    }
  }

  return ctUnion;
}

PolyFillType fill_type( char *f ) {
  int i, j, k;
  char x[][16] = { "EvenOdd", "NonZero", "Positive", "Negative",
                   "evenodd", "nonzero", "positive", "negative",
                   "even-odd", "non-zero", "pos", "neg",
                   "eo", "nz", "p", "n" };
  PolyFillType lookup[16] = { pftEvenOdd, pftNonZero, pftPositive, pftNegative,
                             pftEvenOdd, pftNonZero, pftPositive, pftNegative,
                             pftEvenOdd, pftNonZero, pftPositive, pftNegative,
                             pftEvenOdd, pftNonZero, pftPositive, pftNegative };

  for (i=0; i<16; i++) {
    if (strncmp( x[i], f, 16 )==0) {
      return lookup[i];
    }
  }

  return pftNonZero;

}

// If we get into an edge case of no points, just return 1.
// If the cleaned polygon has a different amount of points
//   than the original, than return (we assume we've been
//   given a cleaned polygon).
// After we clean and simplify the polygon, walk the
//   resulting simplified polygon.  If there is a point
//   that doesn't match the original, return 1.
// Cleaning will remove co-linear points. Assuming co-linear
//   points have been removed, the only way that extra points
//   could be added if there was a self intersection.
//   If we cind an extra point, we know the polygon
//   self intersects.
//
// If there aren't any self intersections, return 0.
//
int self_intersect_test( Path &poly ) {
  Paths opolys;
  int beg_ind, i, j, k, n;

  if (poly.size()<1) { return 1; }

  CleanPolygon(poly);
  SimplifyPolygon( poly, opolys );

  if (poly.size()<1) { return 1; }
  if (opolys.size() != 1) { return 1; }
  if (opolys[0].size() < 1) { return 1; }
  if (opolys[0].size() != poly.size()) { return 1; }

  for (beg_ind=0; beg_ind<opolys[0].size(); beg_ind++) {
    if ( ( opolys[0][beg_ind].X == poly[0].X ) && 
         ( opolys[0][beg_ind].Y == poly[0].Y ) )
      break;
  }


  if (beg_ind == opolys[0].size()) { return 1; }

  n = poly.size();

  for (k=1, i = ((beg_ind+1)%n); i != beg_ind; i = ((i+1)%n), k = ((k+1)%n) ) {
    if ( ( opolys[0][i].X != poly[k].X ) && 
         ( opolys[0][i].Y != poly[k].Y ) )
      return 1;
  }

  return 0;
}

void sort_schedule_delaunay(PolyNode *nod, sort_order_t *sort_order) {
  unsigned long long int X, Y;
  PolyNode *tnod;
  int i, j, k, idx;
  PolyNodes *nodes;


  std::vector<p2t::Point *> pnts;

  nodes = &(nod->Childs);
  if (nodes->size()==0) { return; }

  for (i=0; i<nodes->size(); i++) {

    tnod = nodes->at(i);
    if (tnod->Contour.size()<1) { continue; }
    X = tnod->Contour[0].X;
    Y = tnod->Contour[0].Y;
    for (idx=1; idx<tnod->Contour.size(); idx++) {
      if (tnod->Contour[idx].X < X) {
        X = tnod->Contour[idx].X;
        Y = tnod->Contour[idx].Y;
      } else if ((tnod->Contour[idx].X == X) && (tnod->Contour[idx].Y < Y)) {
        X = tnod->Contour[idx].X;
        Y = tnod->Contour[idx].Y;
      }
    }

    sort_order[i].orig_idx = i;
    sort_order[i].X = X;
    sort_order[i].Y = Y;

    printf("## %lli %lli\n", X, Y);

    pnts.push_back( new p2t::Point( X, Y ) );

  }

  if (pnts.size() > 2) {

    p2t::CDT *cdt = new p2t::CDT(pnts);
    cdt->Triangulate();
    std::vector< p2t::Triangle * > tris = cdt->GetTriangles();

    for (i=0; i<tris.size(); i++) {
      double x, y;
      p2t::Triangle *tri = tris[i];
      x = tri->GetPoint(0)->x;
      y = tri->GetPoint(0)->y;

      printf("%f %f\n", tri->GetPoint(0)->x, tri->GetPoint(0)->y);
      printf("%f %f\n", tri->GetPoint(1)->x, tri->GetPoint(1)->y);
      printf("%f %f\n", tri->GetPoint(2)->x, tri->GetPoint(2)->y);
      printf("\n");
    }
    delete cdt;
  }


  for (i=0; i<pnts.size(); i++) { delete pnts[i]; }

  exit(0);
}

void sort_schedule(PolyNode *nod, sort_order_t *sort_order) {
  unsigned long long int X, Y;
  PolyNode *tnod;
  int i, j, k, idx;
  PolyNodes *nodes;

  nodes = &(nod->Childs);
  if (nodes->size()==0) { return; }

  for (i=0; i<nodes->size(); i++) {

    tnod = nodes->at(i);
    if (tnod->Contour.size()<1) { continue; }
    X = tnod->Contour[0].X;
    Y = tnod->Contour[0].Y;
    for (idx=1; idx<tnod->Contour.size(); idx++) {
      if (tnod->Contour[idx].X < X) {
        X = tnod->Contour[idx].X;
        Y = tnod->Contour[idx].Y;
      } else if ((tnod->Contour[idx].X == X) && (tnod->Contour[idx].Y < Y)) {
        X = tnod->Contour[idx].X;
        Y = tnod->Contour[idx].Y;
      }
    }

    sort_order[i].orig_idx = i;
    sort_order[i].X = X;
    sort_order[i].Y = Y;

  }

  qsort(sort_order, nodes->size(), sizeof(sort_order_t), sort_order_xy_cmp);

}

void walk_poly_tree(PolyNode *nod, int level) {
  unsigned long long int X, Y;
  PolyNode *tnod;
  int i, j, k, idx;
  PolyNodes *nodes;
  sort_order_t *sort_order = NULL;

  if (nod==NULL) { return; }

  if (gSortOrder<0) {
    if (g_verbose_level>0) {
      printf("# %i (IsHole: %i, IsOpen: %i, lvl: %i)\n", (int)(nod->Contour.size()), (int)(nod->IsHole()), (int)(nod->IsOpen()), level );
    }

    if (gPrintReverse) {
      for (i=(nod->Contour.size()-1); i>=0; i--) {
        printf("%lli %lli\n", nod->Contour[i].X, nod->Contour[i].Y);
      }
    } else {
      for (i=0; i<(nod->Contour.size()); i++) {
        printf("%lli %lli\n", nod->Contour[i].X, nod->Contour[i].Y);
      }
    }
    printf("\n");
  }

  nodes = &(nod->Childs);
  if (nodes->size()>0) {

    sort_order = (sort_order_t *)malloc(sizeof(sort_order_t)*nodes->size());
    sort_schedule(nod, sort_order);
    //sort_schedule_delaunay(nod, sort_order);

    /*
    for (i=0; i<nodes->size(); i++) {

      tnod = nodes->at(i);
      if (tnod->Contour.size()<1) { continue; }
      X = tnod->Contour[0].X;
      Y = tnod->Contour[0].Y;
      for (idx=1; idx<tnod->Contour.size(); idx++) {
        if (tnod->Contour[idx].X < X) {
          X = tnod->Contour[idx].X;
          Y = tnod->Contour[idx].Y;
        } else if ((tnod->Contour[idx].X == X) && (tnod->Contour[idx].Y < Y)) {
          X = tnod->Contour[idx].X;
          Y = tnod->Contour[idx].Y;
        }
      }

      sort_order[i].orig_idx = i;
      sort_order[i].X = X;
      sort_order[i].Y = Y;

    }

    qsort(sort_order, nodes->size(), sizeof(sort_order_t), sort_order_xy_cmp);
    */

    for (i=0; i<nodes->size(); i++) {
      //walk_poly_tree(nodes->at(i), level+1);
      walk_poly_tree(nodes->at(sort_order[i].orig_idx), level+1);
    }

    free(sort_order);

  }

  if (gSortOrder>=0) {
    if (g_verbose_level>0) {
      printf("# %i (IsHole: %i, IsOpen: %i, lvl: %i)\n", (int)(nod->Contour.size()), (int)(nod->IsHole()), (int)(nod->IsOpen()), level );
    }
    for (i=0; i<(nod->Contour.size()); i++) {
      printf("%lli %lli\n", nod->Contour[i].X, nod->Contour[i].Y);
    }
    printf("\n");
  }

}

void process_poly_tree(PolyTree &soln_tree) {
  int i;

  //DEBUG
  walk_poly_tree(&soln_tree, 0);
  return;
  //DEBUG

  for (i=0; i<soln_tree.Childs.size(); i++) {
    walk_poly_tree(soln_tree.Childs[i], 0);
  }

}

int main(int argc, char **argv) {
  int i, j, k, idx;
  char ch;
  bool res;
  bool poly_offset_flag = false;

  char ct[][16] = { "Intersection", "Union", "Difference", "Xor" };
  char pt[][16] = { "EvenOdd", "NonZero", "Positive", "Negative" };

  sort_order_t *sort_order = NULL;

  PolyFillType subj_type=pftNonZero, clip_type=pftNonZero;
  ClipType clip_op_type = ctUnion;


  Paths subj_polys, clip_polys, soln;
  PolyTree soln_tree;
  Clipper clip;

  vector<char *> subj_fn;
  vector<char *> clip_fn;

  while ((ch = getopt(argc, argv, "f:s:c:t:S:C:x:vVR:M:E:O:FP:TrBb:")) != -1) {
    switch (ch) {
      case 'f':
        g_function = strdup(optarg);
        break;
      case 's':
        subj_fn.push_back( strdup(optarg) );
        break;
      case 'c':
        clip_fn.push_back( strdup(optarg) );
        break;
      case 't':
        clip_op_type = find_clip_type( optarg );
        break;
      case 'S':
        subj_type = fill_type( optarg );
        break;

      case 'C':
        clip_type = fill_type( optarg );
        break;
      case 'M':
        gOffsetMiterLimit = atof(optarg);
        break;
      case 'E':
        gEps = atof(optarg);
        break;

      case 'B':
        gPrintBoundingBox = true;
        break;
      case 'b':
        gBoundingBoxMargin = atoi(optarg);
        break;
      case 'r':
        gPrintReverse = true;
        break;
      case 'P':
        processSortOrder(optarg);
        break;
      case 'T':
        gPolyTreeFlag = 1;
        break;

      case 'F':
        gReadFloatFlag = 1;
        break;

      case 'R':
        gOffsetRadius = atof(optarg);
        poly_offset_flag = true;
        break;
      case 'x':
        gMulFactor = atoll( optarg );
        if (gMulFactor < 0) gMulFactor = -1*gMulFactor;
        break;
      case 'O':
        gForceOrientation = atoi( optarg );
        break;
      case 'v':
        //g_verbose_flag = true;
        g_verbose_level++;
        break;
      case 'V':
        show_version();
        exit(0);
        break;
      default:
        printf("bad argument\n");
        show_help();
        exit(1);
    }
  }

  if ((subj_fn.size()==0) && (clip_fn.size()==0)) {
    show_help();
    exit(1);
  }

  load_polys( subj_fn, subj_polys );
  load_polys( clip_fn, clip_polys );

  if (g_function && (strncmp(g_function, "self-intersect", strlen("self-intersect"))==0)) {
    for (i=0; i<subj_polys.size(); i++) {
      int b;
      b = self_intersect_test( subj_polys[i] );
      printf("# subj %i self-intersect %i\n", i, b);
    }
    exit(0);
  }

  //if (g_verbose_flag) {
  if (g_verbose_level>1) {
    printf("# clipop: %s, subjtype: %s, cliptype: %s\n", ct[(int)clip_op_type], pt[(int)subj_type], pt[(int)clip_type]);
    printf("# subject (%i)\n", (int)subj_polys.size());
    for (i=0; i<subj_polys.size(); i++) {
      printf("# subj %i\n", i);
      for (j=0; j<subj_polys[i].size(); j++) {
        printf("%lli %lli\n", subj_polys[i][j].X, subj_polys[i][j].Y);
      }
      printf("\n");
    }
    printf("\n");

    printf("# clip (%i)\n", (int)clip_polys.size());
    for (i=0; i<clip_polys.size(); i++) {
      printf("# clip %i\n", i);
      for (j=0; j<clip_polys[i].size(); j++) {
        printf("%lli %lli\n", clip_polys[i][j].X, clip_polys[i][j].Y);
      }
      printf("\n");
    }
  }

  clip.AddPaths( subj_polys, ptSubject, true );
  clip.AddPaths( clip_polys, ptClip, true );

  int count=0;
  PolyNode *nod;

  if (gPolyTreeFlag>0) {
    PolyTree offset_soln_tree;

    res = clip.Execute( clip_op_type, soln_tree, subj_type, clip_type );
    if (!res) { fprintf(stderr, "ERROR\n"); exit(1); }

    if (poly_offset_flag) {
      if (gOffsetRadius > gEps) {

        ClipperOffset co;
        Paths paths;

        nod = soln_tree.GetFirst();

        while (nod) {
          paths.push_back(nod->Contour);
          nod = nod->GetNext();
          count++;
        }

        co.MiterLimit = gOffsetMiterLimit;
        co.AddPaths(paths, jtMiter, etClosedPolygon);

        co.Execute(offset_soln_tree, gMulFactor * gOffsetRadius );

        soln_tree = offset_soln_tree;
      }
    }

    unsigned long long int minX, maxX, minY, maxY;

    process_poly_tree(soln_tree);

    if (gPrintBoundingBox) {
      count=0;
      nod = soln_tree.GetFirst();
      while (nod) {
        if (nod->Contour.size()==0) { continue; }
        if (count==0) {
          minX = nod->Contour[0].X;
          maxX = nod->Contour[0].X;
          minY = nod->Contour[0].Y;
          maxY = nod->Contour[0].Y;
        }
        for (i=0; i<nod->Contour.size(); i++) {
          if (minX > nod->Contour[i].X) { minX = nod->Contour[i].X; }
          if (minY > nod->Contour[i].Y) { minY = nod->Contour[i].Y; }
          if (maxX < nod->Contour[i].X) { maxX = nod->Contour[i].X; }
          if (maxY < nod->Contour[i].Y) { maxY = nod->Contour[i].Y; }
        }
        nod = nod->GetNext();
        count++;
      }

      unsigned long long int marg = gBoundingBoxMargin * gMulFactor;

      printf("%lli %lli\n", minX - marg, minY - marg);
      printf("%lli %lli\n", minX - marg, maxY + marg);
      printf("%lli %lli\n", maxX + marg, maxY + marg);
      printf("%lli %lli\n", maxX + marg, minY - marg);
      printf("\n");

    }




    exit(0);
  }

  res = clip.Execute( clip_op_type, soln, subj_type, clip_type );
  if (!res) { fprintf(stderr, "ERROR\n"); exit(1); }

  //printf("### %i %f\n", poly_offset_flag, gOffsetRadius);

  if (poly_offset_flag) {
    printf("# offsetting... ???\n");
    if (gOffsetRadius > gEps) {

      printf("# offsetting...\n");


      Paths offset_soln;
      ClipperOffset co;
      co.MiterLimit = gOffsetMiterLimit;
      co.AddPaths(soln, jtMiter, etClosedPolygon);
      co.Execute(offset_soln, gMulFactor * gOffsetRadius );
      soln = offset_soln;
    }
  }

  if (gSortOrder>0) {
    sort_order = (sort_order_t *)malloc(sizeof(sort_order_t)*((int)soln.size()));
    for (i=0; i<soln.size(); i++) {
      sort_order[i].orig_idx = i;

      if ((gSortOrder == 1) || (gSortOrder == 2)) {
        sort_order[i].val = ClipperLib::Area(soln[i]);
      } else {
        sort_order[i].val = (double)ClipperLib::Orientation(soln[i]);
      }
    }

    qsort(sort_order, (size_t)soln.size(), sizeof(sort_order_t), sort_order_cmp);

    if (g_verbose_level>1) {
      for (i=0; i<soln.size(); i++) {
        printf("# [%i] orig_idx: %i, val: %f\n", i, sort_order[i].orig_idx, (float)sort_order[i].val);
      }
    }
  }

  if (g_verbose_level>0) { printf("# (%i) sort-order: %i\n", (int)soln.size(), gSortOrder); }
  for (i=0; i<soln.size(); i++) {
    idx = i;
    if (gSortOrder>0) { idx = sort_order[i].orig_idx; }

    if (g_verbose_level>0) { printf("# %i (orig_idx: %i) (orientation: %i)\n", i, idx, (int) ClipperLib::Orientation(soln[idx]) ); }

    for (j=0; j<soln[idx].size(); j++) {
      printf("%lli %lli\n", soln[idx][j].X, soln[idx][j].Y);
    }
    printf("\n");
  }



}
