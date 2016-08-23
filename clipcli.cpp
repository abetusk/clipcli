#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <unistd.h>
#include <string.h>

#include <errno.h>

#include <vector>

#include "clipper.hpp"

using namespace std;
using namespace ClipperLib;

bool g_verbose_flag=false;
char gVersion[]="0.0.1";

char *g_function=NULL;

long long int gMulFactor=1;
double gOffsetMiterLimit = 3.0;
double gOffsetRadius = 0.0;
double gEps = 0.000001;


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
  printf("  [-v]                Verbose.\n");
  printf("  [-V]                Show version.\n");
  printf("\n");

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

    k = sscanf(buf, "%lli %lli", &X, &Y);
    if (k!=2) { fprintf(stderr, "invalid read on line %i, exiting\n", line_no); exit(2); }

    cur_clip_path.push_back( ClipperLib::IntPoint(X*gMulFactor, Y*gMulFactor) );
  }

  if (cur_clip_path.size() > 0) {
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

int main(int argc, char **argv) {
  int i, j, k;
  char ch;
  bool res;
  bool poly_offset_flag = false;

  char ct[][16] = { "Intersection", "Union", "Difference", "Xor" };
  char pt[][16] = { "EvenOdd", "NonZero", "Positive", "Negative" };

  PolyFillType subj_type=pftNonZero, clip_type=pftNonZero;
  ClipType clip_op_type = ctUnion;

  Paths subj_polys, clip_polys, soln;
  Clipper clip;

  vector<char *> subj_fn;
  vector<char *> clip_fn;

  while ((ch = getopt(argc, argv, "f:s:c:t:S:C:x:vVR:M:E:")) != -1) {
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

      case 'R':
        gOffsetRadius = atof(optarg);
        poly_offset_flag = true;
        break;
      case 'x':
        gMulFactor = atoll( optarg );
        if (gMulFactor < 0) gMulFactor = -1*gMulFactor;
        break;
      case 'v':
        g_verbose_flag = true;
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

  if (g_verbose_flag) {
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

  res = clip.Execute( clip_op_type, soln, subj_type, clip_type );
  if (!res) { fprintf(stderr, "ERROR\n"); exit(1); }

  printf("### %i %f\n", poly_offset_flag, gOffsetRadius);

  if (poly_offset_flag) {
    if (gOffsetRadius > gEps) {



      Paths offset_soln;
      ClipperOffset co;
      co.MiterLimit = gOffsetMiterLimit;
      co.AddPaths(soln, jtMiter, etClosedPolygon);
      co.Execute(offset_soln, gMulFactor * gOffsetRadius );
      soln = offset_soln;
    }
  }

  printf("# (%i)\n", (int)soln.size());
  for (i=0; i<soln.size(); i++) {
    printf("# %i\n", i);
    for (j=0; j<soln[i].size(); j++) {
      printf("%lli %lli\n", soln[i][j].X, soln[i][j].Y);
    }
    printf("\n");
  }



}
