# clipcli

A command line tool to play around with [clipperlib](http://www.angusj.com/delphi/clipper.php).

usage:

```
clipcli [-t clip-type] [-S subj-fill] [-C clip-fill] [-v] [-V] [-s subj] [-s subj] ... [-c clip] [-c clip] ...

  [-t clip-type]      Clipping operation.  One of (Intersection, Union, Difference, Xor).  Union default.
  [-S subj-fill]      Subject fill type.  One of (EvenOdd, NonZero, Positive, Negative).  NonZero default.
  [-C clip-fill]      Clip fill type.  One of (EvenOdd, NonZero, Positive, Negative).  NonZero default.
  [-s subj]           File containg subject polygon.
  [-c clip]           File containg clip polygon.
  [-x mul]            Apply multiplication factor to input polygons.
  [-R radius]         Polygon offset radius.
  [-M miter_limit]    Miter limit (default 3.000000).
  [-E epsilon]        Epsilon (default 0.000001).
  [-O orientation]    force orientation of all polygons (>0 cw, <0 ccw, 0 (default) input orientation).
  [-F]                read floating point input instead of long long integer
  [-v]                Increase verbose level (can be specified multiple times).
  [-P sortorder]      Polygon output order.  One of (area-inc|area-dec|cw-inc|cw-dec|pre|post|undef|infix|closest).
  [-B]                Print bounding box (cw)
  [-b margin]         Add margin to bounding box
  [-r]                print reverse order polygon contours
  [-T]                Output in tree order
  [-V]                Show version.
```

# Example

```bash
$ cat testdata/colinearbox.gp ; echo ; ./clipcli -s testdata/colinearbox.gp 
0 0
50 0
100 0
100 50
100 100
50 100
0 100
0 50

100 100
0 100
0 0
100 0

```

# Credits

* [Angus Johnson's ClipperLib](http://www.angusj.com/delphi/clipper.php)
* [Bl4ckb0ne's delaunay-triangulation](https://github.com/Bl4ckb0ne/delaunay-triangulation)

# License

Unless otherwise specified, the license is AGPLv3.

All libraries have been chosen to have a GPL compatible license.
See individual files for license details.
