/*
 *  clipcli is a polygon clipping tool that acts as a front end to Angus Johnson's clipperlib library
 *  Copyright (C) 2019 abetusk
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Affero General Public License as
 *  published by the Free Software Foundation, either version 3 of the
 *  License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Affero General Public License for more details.
 *
 *  You should have received a copy of the GNU Affero General Public License
 *  along with this program.  If not, see <https://www.gnu.org/licenses/>.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <unistd.h>
#include <string.h>

#include <errno.h>

#include <vector>
#include <map>

#include "clipper.hpp"

#include "vector2.h"
#include "triangle.h"
#include "delaunay.h"

char gVersion[]="0.1.1";

int g_verbose_level = 0;
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

bool gConvexHull = false;

using namespace std;
using namespace ClipperLib;

typedef signed long long cInt;
typedef unsigned long long cUInt;

typedef Vector2<double> Vec2f;

cInt dtocint( double d ) {
  if (d < 0.0) return (unsigned long long)(d-0.5);
  return (signed long long)(d+0.5);
}


class iPnt {
  public:
  long long int X, Y;
  iPnt() { X=0; Y=0; }
  iPnt(long long int a, long long int b) { X = a; Y = b; }
} ;

struct iPnt_cmp {
  bool operator() (const iPnt &lhs, const iPnt &rhs) const {
    if ( lhs.X == rhs.X ) return lhs.Y < rhs.Y ;
    return lhs.X < rhs.X;
  }
};

class PointInfo {
  public:
  int pointInd, polygonInd;
  int dir;
  bool visited, is_outer_boundary, is_clockwise;
} ;

typedef std::map< iPnt , PointInfo, iPnt_cmp > PointInfoMap;
typedef std::map< iPnt, std::vector< iPnt >, iPnt_cmp > PointAdjacency;


class iEdge {
  public:
    iPnt v[2];
};

struct iEdge_cmp {
  bool operator() (const iEdge &lhs, const iEdge &rhs) const {
    if ( ( lhs.v[0].X == rhs.v[0].X ) &&
         ( lhs.v[0].Y == rhs.v[0].Y ) ) {
      if ( lhs.v[1].X == rhs.v[1].X ) return lhs.v[1].Y < rhs.v[1].Y;
      return lhs.v[1].X < rhs.v[1].X;
    }

    if ( lhs.v[0].X == rhs.v[0].X ) return lhs.v[0].Y < rhs.v[0].Y ;
    return lhs.v[0].X < rhs.v[0].X;
  }
};

typedef std::map< iEdge, int, iEdge_cmp > EdgeMap;


void addEdgesFromTri( PointAdjacency &edge, Triangle *tri, EdgeMap &edgeMap ){
  iEdge uv;
  int a,b,i, perm[][2] = { {0,1}, {1,0}, {0,2}, {2,0}, {1,2}, {2,1} };

  uv.v[0].X = dtocint( tri->p1.x ); uv.v[0].Y = dtocint( tri->p1.y );
  uv.v[1].X = dtocint( tri->p2.x ); uv.v[1].Y = dtocint( tri->p2.y );

  if ( edgeMap.find( uv ) == edgeMap.end() ) {

    //printf("ADDING (%lli %lli) (%lli %lli)\n", uv.v[0].X, uv.v[0].Y, uv.v[1].X, uv.v[1].Y);

    edge[ uv.v[0] ].push_back( uv.v[1] );
    edgeMap[ uv ] = 1;
  }

  uv.v[0].X = dtocint( tri->p2.x ); uv.v[0].Y = dtocint( tri->p2.y );
  uv.v[1].X = dtocint( tri->p1.x ); uv.v[1].Y = dtocint( tri->p1.y );

  if ( edgeMap.find( uv ) == edgeMap.end() ) {

    //printf("ADDING (%lli %lli) (%lli %lli)\n", uv.v[0].X, uv.v[0].Y, uv.v[1].X, uv.v[1].Y);

    edge[ uv.v[0] ].push_back( uv.v[1] );
    edgeMap[ uv ] = 1;
  }

  //--

  uv.v[0].X = dtocint( tri->p1.x ); uv.v[0].Y = dtocint( tri->p1.y );
  uv.v[1].X = dtocint( tri->p3.x ); uv.v[1].Y = dtocint( tri->p3.y );

  if ( edgeMap.find( uv ) == edgeMap.end() ) {

    //printf("ADDING (%lli %lli) (%lli %lli)\n", uv.v[0].X, uv.v[0].Y, uv.v[1].X, uv.v[1].Y);

    edge[ uv.v[0] ].push_back( uv.v[1] );
    edgeMap[ uv ] = 1;
  }

  uv.v[0].X = dtocint( tri->p3.x ); uv.v[0].Y = dtocint( tri->p3.y );
  uv.v[1].X = dtocint( tri->p1.x ); uv.v[1].Y = dtocint( tri->p1.y );

  if ( edgeMap.find( uv ) == edgeMap.end() ) {

    //printf("ADDING (%lli %lli) (%lli %lli)\n", uv.v[0].X, uv.v[0].Y, uv.v[1].X, uv.v[1].Y);

    edge[ uv.v[0] ].push_back( uv.v[1] );
    edgeMap[ uv ] = 1;
  }

  //--

  uv.v[0].X = dtocint( tri->p2.x ); uv.v[0].Y = dtocint( tri->p2.y );
  uv.v[1].X = dtocint( tri->p3.x ); uv.v[1].Y = dtocint( tri->p3.y );

  if ( edgeMap.find( uv ) == edgeMap.end() ) {

    //printf("ADDING (%lli %lli) (%lli %lli)\n", uv.v[0].X, uv.v[0].Y, uv.v[1].X, uv.v[1].Y);

    edge[ uv.v[0] ].push_back( uv.v[1] );
    edgeMap[ uv ] = 1;
  }

  uv.v[0].X = dtocint( tri->p3.x ); uv.v[0].Y = dtocint( tri->p3.y );
  uv.v[1].X = dtocint( tri->p2.x ); uv.v[1].Y = dtocint( tri->p2.y );

  if ( edgeMap.find( uv ) == edgeMap.end() ) {

    //printf("ADDING (%lli %lli) (%lli %lli)\n", uv.v[0].X, uv.v[0].Y, uv.v[1].X, uv.v[1].Y);

    edge[ uv.v[0] ].push_back( uv.v[1] );
    edgeMap[ uv ] = 1;
  }

}




#define SORT_DELAUNAY 16

void show_version(void) {
  printf("%s\n", gVersion);
}


void show_help(void) {
  printf("version: "); show_version(); printf("\n");
  printf("usage:\n\n");
  printf("    clipcli [-t clip-type] [-S subj-fill] [-C clip-fill] [-v] [-V] [-s subj] [-s subj] ... [-c clip] [-c clip] ...\n\n");
  printf("  [-t clip-type]      Clipping operation.  One of (Intersection, Union, Difference, Xor).  Union default.\n");
  printf("  [-S subj-fill]      Subject fill type.  One of (EvenOdd, NonZero, Positive, Negative).  NonZero default.\n");
  printf("  [-C clip-fill]      Clip fill type.  One of (EvenOdd, NonZero, Positive, Negative).  NonZero default.\n");
  printf("  [-s subj]           File containing subject polygon.\n");
  printf("  [-c clip]           File containing clip polygon.\n");
  printf("  [-x mul]            Apply multiplication factor to input polygons.\n");
  printf("  [-R radius]         Polygon offset radius (negative for infill).\n");
  printf("  [-M miter_limit]    Miter limit (default %f).\n", gOffsetMiterLimit);
  printf("  [-E epsilon]        Epsilon (default %f).\n", gEps);
  printf("  [-O orientation]    force orientation of all polygons (>0 cw, <0 ccw, 0 (default) input orientation).\n");
  printf("  [-F]                read floating point input instead of long long integer\n");
  printf("  [-v]                Increase verbose level (can be specified multiple times).\n");
  printf("  [-P sortorder]      Polygon output order.  One of (area-inc|area-dec|cw-inc|cw-dec|pre|post|undef|infix|closest).\n");
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
  } else if (strcmp(so, "closest")==0) {
    gSortOrder = SORT_DELAUNAY;
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

int sort_order_cmp_2(const void *x, const void *y) {
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


void load_poly_set_orientation(Path &cur_clip_path) {

  if (cur_clip_path.size() > 0) {

    //printf("##?? gForceOrientation %i\n", gForceOrientation);

    if (gForceOrientation!=0) {

      if (gForceOrientation>0) {
        if (!ClipperLib::Orientation(cur_clip_path)) {

          //printf("##>> gForceOrientation > 0 (%i), reversing...\n", gForceOrientation);

          ClipperLib::ReversePath(cur_clip_path);
        }
      } else {
        if (ClipperLib::Orientation(cur_clip_path)) {

          //printf("##<< gForceOrientation < 0 (%i), reversing...\n", gForceOrientation);

          ClipperLib::ReversePath(cur_clip_path);
        }
      }

    }

  }

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
      load_poly_set_orientation(cur_clip_path);
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
    load_poly_set_orientation(cur_clip_path);
    p.push_back( cur_clip_path );
  }

  /*
  printf("##?? ... gForceOrientation %i\n", gForceOrientation);

  if (cur_clip_path.size() > 0) {

    printf("##?? gForceOrientation %i\n", gForceOrientation);

    if (gForceOrientation!=0) {

      if (gForceOrientation>0) {
        if (!ClipperLib::Orientation(cur_clip_path)) {

          printf("##>> gForceOrientation > 0 (%i), reversing...\n", gForceOrientation);

          ClipperLib::ReversePath(cur_clip_path);
        }
      } else {
        if (ClipperLib::Orientation(cur_clip_path)) {

          printf("##<< gForceOrientation < 0 (%i), reversing...\n", gForceOrientation);

          ClipperLib::ReversePath(cur_clip_path);
        }
      }

    }

    p.push_back(cur_clip_path);

  }
  */

}


void load_polys( vector<char *> &fns, Paths &polys) {
  FILE *fp;
  int i, j, k;

  for (i=0; i<fns.size(); i++) {
    fp = fopen(fns[i], "r");
    if (fp==NULL) {
      perror(fns[i]);
      exit(-1);
    }
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

int walk_recur_sort_order_r(sort_order_t *sort_order, int pos, int cur_order, std::vector< std::vector<int> > &edge) {
  int i, j, k, orig_order;
  iPnt p, q;
  std::vector<int> v;

  if (sort_order[pos].val>=0) { return 0; }

  orig_order = cur_order;

  sort_order[pos].val = cur_order;
  cur_order++;

  v = edge[pos];
  for (i=0; i<v.size(); i++) {
    cur_order += walk_recur_sort_order_r(sort_order, v[i], cur_order, edge);
  }

  return cur_order - orig_order;
}

void walk_recur_sort_order(sort_order_t *sort_order, std::vector< std::vector<int> > &edge) {
  int i, j, cur_order, n;
  n = (int)(edge.size());
  cur_order=0;
  for (i=0; i<n; i++) {
    if (sort_order[i].val>=0) { continue; }
    cur_order += walk_recur_sort_order_r(sort_order, i, cur_order, edge);
  }

}


void sort_schedule_delaunay(PolyNode *nod, sort_order_t *sort_order) {
  unsigned long long int X, Y;
  PolyNode *tnod;
  int i, j, k, idx, tot=0;
  PolyNodes *nodes;

  double x, y;
  int local_debug = 0;

  std::vector<Vec2f> pnts;

  nodes = &(nod->Childs);
  if (nodes->size()==0) { return; }

  if (local_debug) {
    printf("\n### sort_schedule_delaunay (%i)\n", (int)(nodes->size()));
  }



  for (i=0; i<nodes->size(); i++) {

    sort_order[i].orig_idx = i;
    sort_order[i].X = 0;
    sort_order[i].Y = 0;
    sort_order[i].val = -2;

    tnod = nodes->at(i);
    if (tnod->Contour.size()<1) {

      if (local_debug) {
        printf("### skipping i %i (contour %i)\n", i, (int)(tnod->Contour.size()));
      }

      continue;
    }
    X = tnod->Contour[0].X;
    Y = tnod->Contour[0].Y;

    sort_order[i].X = X;
    sort_order[i].Y = Y;
    sort_order[i].val = -1;

    if (local_debug) {
      printf("## %lli %lli (%f)\n", X, Y, sort_order[i].val);
    }

    pnts.push_back( Vec2f((double)X, (double)Y) );

  }

  if (nodes->size() < 3) { return; }


  if (local_debug) {
    printf("\n\n### pnts %i\n", (int)(pnts.size()));
    for (i=0; i<pnts.size(); i++) {
      printf("#### %f %f\n", pnts[i].x, pnts[i].y);
    }
    printf("\n\n");
  }

  Delaunay triangulation;
  std::vector<Triangle> tris = triangulation.triangulate(pnts);
  std::vector<Edge> edges = triangulation.getEdges();

  if (local_debug) {
    for (i=0; i<edges.size(); i++) {
      printf("##_ %f %f\n##_ %f %f\n\n",
          edges[i].p1.x, edges[i].p1.y,
          edges[i].p2.x, edges[i].p2.y);
    }
  }

  PointAdjacency pointAdj;
  EdgeMap edgeMap;

  for (i=0; i<tris.size(); i++) {
    addEdgesFromTri(pointAdj, &(tris[i]), edgeMap);
  }

  std::vector< std::vector<int> > adj;
  std::vector< int > vv;
  std::vector< iPnt > pvec;
  PointInfo pinfo;
  PointInfoMap pim;
  iPnt p;

  for (i=0; i<nodes->size(); i++) {
    if (sort_order[i].val == -2) { continue; }

    p.X = sort_order[i].X;
    p.Y = sort_order[i].Y;

    if (pointAdj.find(p) == pointAdj.end()) {
      printf("!! ERROR, COULD NOT FIND %lli %lli INT pointAdj\n", p.X, p.Y);
      continue;
    }


    pinfo.pointInd = i;
    pim[p] = pinfo;

    if (local_debug) {
      printf("### adding %lli %lli idx %i\n", p.X, p.Y, pinfo.pointInd);
    }

  }


  for (i=0; i<nodes->size(); i++) {
    if (sort_order[i].val == -2) { continue; }

    p.X = sort_order[i].X;
    p.Y = sort_order[i].Y;

    if (pointAdj.find(p) == pointAdj.end()) {
      printf("ERROR, COULD NOT FIND %lli %lli INT pointAdj\n", p.X, p.Y);
      continue;
    }

    vv.clear();
    pvec = pointAdj[p];
    for (j=0; j<pvec.size(); j++) {
      p.X = pvec[j].X;
      p.Y = pvec[j].Y;
      pinfo = pim[p];

      vv.push_back(pinfo.pointInd);
    }

    adj.push_back(vv);

  }

  if (local_debug) {
    for (i=0; i<adj.size(); i++) {
      vv = adj[i];
      for (j=0; j<vv.size(); j++) {
        printf("#_ %lli %lli\n#_ %lli %lli\n#_\n\n",
            sort_order[i].X, sort_order[i].Y,
            sort_order[vv[j]].X, sort_order[vv[j]].Y);

      }
    }
  }

  walk_recur_sort_order(sort_order, adj);

  int n = nodes->size();
  if (local_debug) {
    printf("\n\n");
    for (i=0; i<n; i++) {
      printf("#- [%i] %lli %lli order: %f\n", i, sort_order[i].X, sort_order[i].Y, sort_order[i].val);
    }
    printf("\n\n");
  }

  qsort(sort_order, n, sizeof(sort_order_t), sort_order_cmp);

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

    if (gSortOrder == SORT_DELAUNAY) {
      sort_schedule_delaunay(nod, sort_order);
    } else {
      sort_schedule(nod, sort_order);
    }

    for (i=0; i<nodes->size(); i++) {
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

  walk_poly_tree(&soln_tree, 0);

  //for (i=0; i<soln_tree.Childs.size(); i++) { walk_poly_tree(soln_tree.Childs[i], 0); }

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

  Path hull_points, hull_path;

  while ((ch = getopt(argc, argv, "f:s:c:t:S:C:x:vVR:M:E:O:FP:TrBb:H")) != -1) {
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

      case 'H':
        gConvexHull = true;
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

  if (gConvexHull) {
    for (i=0; i<subj_polys.size(); i++) {
      for (j=0; j<subj_polys[i].size(); j++) {
        hull_points.push_back( subj_polys[i][j] );
      }
    }

    ConvexHull(hull_points, hull_path);

    for (i=0; i<hull_path.size(); i++) {
      printf("%lli %lli\n", hull_path[i].X, hull_path[i].Y);
    }

    exit(0);
  }

  load_polys( clip_fn, clip_polys );

  if (g_function && (strncmp(g_function, "self-intersect", strlen("self-intersect"))==0)) {
    for (i=0; i<subj_polys.size(); i++) {
      int b;
      b = self_intersect_test( subj_polys[i] );
      printf("# subj %i self-intersect %i\n", i, b);
    }
    exit(0);
  }

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
      if (fabs(gOffsetRadius) > gEps) {

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

    if (fabs(gOffsetRadius) > gEps) {
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
