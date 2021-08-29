
'''
Script to generate an uniform k-grid in the Brillouin zone (in crystal coordinates: K_POINTS (crystal) )
including or not including the periodic images of the Gamma point.

It generates a uniform k-grid inside the unit cube [0,1]**3 in units of the reciprocal lattice vectors bg with the possibility to include the corners of the cube. 
To be used in the input of pw.x for the nscf calculation needed to run simple.x.
'''


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, prog='gen_kgrid',description="""
This script generates an uniform k-grid in the Brillouin zone (in crystal coordinates) 
w/ and w/o the periodic images of the Gamma point. 
To be used in the input of pw.x for the nscf calculation needed for simple.x.

Example:
$ python gen_kgrid.py 12 12 12 --images    --> it generates a 12x12x12 grid including the periodic images of Gamma 
""")
    parser.add_argument("kx", type=int, help="k-grid (x)")
    parser.add_argument("ky", type=int, help="k-grid (y)")
    parser.add_argument("kz", type=int, help="k-grid (z)")
    parser.add_argument("--images", help="include periodic images of Gamma point", action="store_true")
    args = parser.parse_args()

    nkpoints = [args.kx, args.ky, args.kz]

    print 'K_POINTS (crystal)'

    if not args.images:
        print nkpoints[0]*nkpoints[1]*nkpoints[2]
    else:
        print nkpoints[0]*nkpoints[1]*nkpoints[2] + 7

    for i in xrange(nkpoints[0]):
        for j in xrange(nkpoints[1]):
            for k in xrange(nkpoints[2]):
                print float(i)/nkpoints[0], float(j)/nkpoints[1], float(k)/nkpoints[2], 1.0

    if args.images:
        print 0.0, 0.0, 1.0, 1.0
        print 0.0, 1.0, 0.0, 1.0
        print 1.0, 0.0, 0.0, 1.0
        print 1.0, 1.0, 0.0, 1.0
        print 1.0, 0.0, 1.0, 1.0
        print 0.0, 1.0, 1.0, 1.0
        print 1.0, 1.0, 1.0, 1.0


