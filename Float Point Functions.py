def padding (src, left=0, top=0, right=0, bottom=0):
    core   = vs.get_core ()
    w      = src.width
    h      = src.height
    clip   = core.fmtc.resample (src, w+left+right, h+top+bottom, -left, -top, w+left+right, h+top+bottom, kernel="point", fulls=True, fulld=True)
    return clip

def Median (src):
    core     = vs.get_core ()
    def kernel (n, f):
        fout = f.copy ()
    
        for p in range (fout.format.num_planes):
            plane        = f.get_read_array (p)
            dst_plane    = fout.get_write_array (p)
    
            plane_height = len (plane)
            plane_width  = len (plane[0])

            members      = [plane[0][0]] * 9 
       
            for y in range (1, plane_height - 1):
                for x in range (1, plane_width - 1):
                    members[0] = plane[y-1][x-1]
                    members[1] = plane[y][x-1]
                    members[2] = plane[y+1][x-1]
                    members[3] = plane[y-1][x]
                    members[4] = plane[y][x]
                    members[5] = plane[y+1][x]
                    members[6] = plane[y-1][x+1]
                    members[7] = plane[y][x+1]
                    members[8] = plane[y+1][x+1]
                    members.sort ()
                    dst_plane[y][x] = members[4]
        return fout
    clip = padding (src, 1, 1, 1, 1)
    clip = core.std.ModifyFrame (clip=clip, clips=clip, selector=kernel)
    clip = core.std.CropRel (clip, 1, 1, 1, 1)
    return clip

def Deflate (src):
    core     = vs.get_core ()
    def mean (n, f):
        fout = f.copy ()
    
        for p in range (fout.format.num_planes):
            plane        = f.get_read_array (p)
            dst_plane    = fout.get_write_array (p)
    
            plane_height = len (plane)
            plane_width  = len (plane[0])
       
            for y in range (1, plane_height - 1):
                for x in range (1, plane_width - 1):
                    dst_plane[y][x] = (plane[y+1][x+1] + plane[y][x+1] + plane[y-1][x+1] + plane[y+1][x] + plane[y-1][x] + plane[y+1][x-1] + plane[y][x-1] + plane[y-1][x-1]) / 8
        return fout
    avrg = padding (src, 1, 1, 1, 1)
    avrg = core.std.ModifyFrame (clip=avrg, clips=avrg, selector=mean)
    avrg = core.std.CropRel (avrg, 1, 1, 1, 1)
    clip = core.std.Expr ([src, avrg], ["y x < y x ?"])
    return clip

def Inflate (src):
    core     = vs.get_core ()
    def mean (n, f):
        fout = f.copy ()
    
        for p in range (fout.format.num_planes):
            plane        = f.get_read_array (p)
            dst_plane    = fout.get_write_array (p)
    
            plane_height = len (plane)
            plane_width  = len (plane[0])
       
            for y in range (1, plane_height - 1):
                for x in range (1, plane_width - 1):
                    dst_plane[y][x] = (plane[y+1][x+1] + plane[y][x+1] + plane[y-1][x+1] + plane[y+1][x] + plane[y-1][x] + plane[y+1][x-1] + plane[y][x-1] + plane[y-1][x-1]) / 8
        return fout
    avrg = padding (src, 1, 1, 1, 1)
    avrg = core.std.ModifyFrame (clip=avrg, clips=avrg, selector=mean)
    avrg = core.std.CropRel (avrg, 1, 1, 1, 1)
    clip = core.std.Expr ([src, avrg], ["y x > y x ?"])
    return clip

def Minimum (src):
    core     = vs.get_core ()
    def kernel (n, f):
        fout = f.copy ()
    
        for p in range (fout.format.num_planes):
            plane        = f.get_read_array (p)
            dst_plane    = fout.get_write_array (p)
    
            plane_height = len (plane)
            plane_width  = len (plane[0])

            members      = [plane[0][0]] * 8 
       
            for y in range (1, plane_height - 1):
                for x in range (1, plane_width - 1):
                    members[0] = plane[y-1][x-1]
                    members[1] = plane[y][x-1]
                    members[2] = plane[y+1][x-1]
                    members[3] = plane[y-1][x]
                    members[4] = plane[y+1][x]
                    members[5] = plane[y-1][x+1]
                    members[6] = plane[y][x+1]
                    members[7] = plane[y+1][x+1]
                    members.sort ()
                    dst_plane[y][x] = members[0]
        return fout
    clip = padding (src, 1, 1, 1, 1)
    clip = core.std.ModifyFrame (clip=clip, clips=clip, selector=kernel)
    clip = core.std.CropRel (clip, 1, 1, 1, 1)
    return clip

def Maximum (src):
    core     = vs.get_core ()
    def kernel (n, f):
        fout = f.copy ()
    
        for p in range (fout.format.num_planes):
            plane        = f.get_read_array (p)
            dst_plane    = fout.get_write_array (p)
    
            plane_height = len (plane)
            plane_width  = len (plane[0])

            members      = [plane[0][0]] * 8 
       
            for y in range (1, plane_height - 1):
                for x in range (1, plane_width - 1):
                    members[0] = plane[y-1][x-1]
                    members[1] = plane[y][x-1]
                    members[2] = plane[y+1][x-1]
                    members[3] = plane[y-1][x]
                    members[4] = plane[y+1][x]
                    members[5] = plane[y-1][x+1]
                    members[6] = plane[y][x+1]
                    members[7] = plane[y+1][x+1]
                    members.sort ()
                    dst_plane[y][x] = members[7]
        return fout
    clip = padding (src, 1, 1, 1, 1)
    clip = core.std.ModifyFrame (clip=clip, clips=clip, selector=kernel)
    clip = core.std.CropRel (clip, 1, 1, 1, 1)
    return clip

def RemoveGrain (src, mode=11):
    core     = vs.get_core ()
    def rg0 (n, f):
        fout = f.copy ()
    
        for p in range (fout.format.num_planes):
            plane        = f.get_read_array (p)
            dst_plane    = fout.get_write_array (p)
    
            plane_height = len (plane)
            plane_width  = len (plane[0])
       
            for y in range (1, plane_height - 1):
                for x in range (1, plane_width - 1):
                    dst_plane[y][x] = plane[y][x]
        return fout

    def rg1 (n, f):
        fout = f.copy ()
    
        for p in range (fout.format.num_planes):
            plane        = f.get_read_array (p)
            dst_plane    = fout.get_write_array (p)
    
            plane_height = len (plane)
            plane_width  = len (plane[0])

            members      = [plane[0][0]] * 8 
       
            for y in range (1, plane_height - 1):
                for x in range (1, plane_width - 1):
                    members[0] = plane[y-1][x-1]
                    members[1] = plane[y][x-1]
                    members[2] = plane[y+1][x-1]
                    members[3] = plane[y-1][x]
                    members[4] = plane[y+1][x]
                    members[5] = plane[y-1][x+1]
                    members[6] = plane[y][x+1]
                    members[7] = plane[y+1][x+1]
                    members.sort ()
                    minnbr     = members[0]
                    maxnbr     = members[7]
                    dst_plane[y][x] = max (min (plane[y][x], maxnbr), minnbr)
        return fout

    def rg2 (n, f):
        fout = f.copy ()
    
        for p in range (fout.format.num_planes):
            plane        = f.get_read_array (p)
            dst_plane    = fout.get_write_array (p)
    
            plane_height = len (plane)
            plane_width  = len (plane[0])

            members      = [plane[0][0]] * 8 
       
            for y in range (1, plane_height - 1):
                for x in range (1, plane_width - 1):
                    members[0] = plane[y-1][x-1]
                    members[1] = plane[y][x-1]
                    members[2] = plane[y+1][x-1]
                    members[3] = plane[y-1][x]
                    members[4] = plane[y+1][x]
                    members[5] = plane[y-1][x+1]
                    members[6] = plane[y][x+1]
                    members[7] = plane[y+1][x+1]
                    members.sort ()
                    minnbr     = members[1]
                    maxnbr     = members[6]
                    dst_plane[y][x] = max (min (plane[y][x], maxnbr), minnbr)
        return fout

    def rg3 (n, f):
        fout = f.copy ()
    
        for p in range (fout.format.num_planes):
            plane        = f.get_read_array (p)
            dst_plane    = fout.get_write_array (p)
    
            plane_height = len (plane)
            plane_width  = len (plane[0])

            members      = [plane[0][0]] * 8 
       
            for y in range (1, plane_height - 1):
                for x in range (1, plane_width - 1):
                    members[0] = plane[y-1][x-1]
                    members[1] = plane[y][x-1]
                    members[2] = plane[y+1][x-1]
                    members[3] = plane[y-1][x]
                    members[4] = plane[y+1][x]
                    members[5] = plane[y-1][x+1]
                    members[6] = plane[y][x+1]
                    members[7] = plane[y+1][x+1]
                    members.sort ()
                    minnbr     = members[2]
                    maxnbr     = members[5]
                    dst_plane[y][x] = max (min (plane[y][x], maxnbr), minnbr)
        return fout

    def rg4 (n, f):
        fout = f.copy ()
    
        for p in range (fout.format.num_planes):
            plane        = f.get_read_array (p)
            dst_plane    = fout.get_write_array (p)
    
            plane_height = len (plane)
            plane_width  = len (plane[0])

            members      = [plane[0][0]] * 8 
       
            for y in range (1, plane_height - 1):
                for x in range (1, plane_width - 1):
                    members[0] = plane[y-1][x-1]
                    members[1] = plane[y][x-1]
                    members[2] = plane[y+1][x-1]
                    members[3] = plane[y-1][x]
                    members[4] = plane[y+1][x]
                    members[5] = plane[y-1][x+1]
                    members[6] = plane[y][x+1]
                    members[7] = plane[y+1][x+1]
                    members.sort ()
                    minnbr     = members[3]
                    maxnbr     = members[4]
                    dst_plane[y][x] = max (min (plane[y][x], maxnbr), minnbr)
        return fout

    def rg5 (n, f):
        fout = f.copy ()
    
        for p in range (fout.format.num_planes):
            plane        = f.get_read_array (p)
            dst_plane    = fout.get_write_array (p)
    
            plane_height = len (plane)
            plane_width  = len (plane[0])
     
            for y in range (1, plane_height - 1):
                for x in range (1, plane_width - 1):
                    max1 = max (plane[y-1][x-1], plane[y+1][x+1])
                    min1 = min (plane[y-1][x-1], plane[y+1][x+1])
                    max2 = max (plane[y-1][x], plane[y+1][x])
                    min2 = min (plane[y-1][x], plane[y+1][x])
                    max3 = max (plane[y-1][x+1], plane[y+1][x-1])
                    min3 = min (plane[y-1][x+1], plane[y+1][x-1])
                    max4 = max (plane[y][x-1], plane[y][x+1])
                    min4 = min (plane[y][x-1], plane[y][x+1])
                    r4   = max (min (plane[y][x], max4), min4)
                    r3   = max (min (plane[y][x], max3), min3)
                    r2   = max (min (plane[y][x], max2), min2)
                    r1   = max (min (plane[y][x], max1), min1)
                    c1   = abs (plane[y][x] - r1)
                    c2   = abs (plane[y][x] - r2)
                    c3   = abs (plane[y][x] - r3)
                    c4   = abs (plane[y][x] - r4)
                    mind = min (min (min (c1, c2), c3), c4)
                    if mind == c4:
                       clip = r4
                    elif mind == c3:
                       clip = r3
                    elif mind == c2:
                       clip = r2
                    else:
                       clip = r1
                    dst_plane[y][x] = clip
        return fout

    def rg6 (n, f):
        fout = f.copy ()
    
        for p in range (fout.format.num_planes):
            plane        = f.get_read_array (p)
            dst_plane    = fout.get_write_array (p)
    
            plane_height = len (plane)
            plane_width  = len (plane[0])
     
            for y in range (1, plane_height - 1):
                for x in range (1, plane_width - 1):
                    max1 = max (plane[y-1][x-1], plane[y+1][x+1])
                    min1 = min (plane[y-1][x-1], plane[y+1][x+1])
                    max2 = max (plane[y-1][x], plane[y+1][x])
                    min2 = min (plane[y-1][x], plane[y+1][x])
                    max3 = max (plane[y-1][x+1], plane[y+1][x-1])
                    min3 = min (plane[y-1][x+1], plane[y+1][x-1])
                    max4 = max (plane[y][x-1], plane[y][x+1])
                    min4 = min (plane[y][x-1], plane[y][x+1])
                    d1   = max1 - min1
                    d2   = max2 - min2
                    d3   = max3 - min3
                    d4   = max4 - min4
                    r4   = max (min (plane[y][x], max4), min4)
                    r3   = max (min (plane[y][x], max3), min3)
                    r2   = max (min (plane[y][x], max2), min2)
                    r1   = max (min (plane[y][x], max1), min1)
                    c1   = max (min ((abs (plane[y][x] - r1) * 2 + d1), 1.0), 0.0)
                    c2   = max (min ((abs (plane[y][x] - r2) * 2 + d2), 1.0), 0.0)
                    c3   = max (min ((abs (plane[y][x] - r3) * 2 + d3), 1.0), 0.0)
                    c4   = max (min ((abs (plane[y][x] - r4) * 2 + d4), 1.0), 0.0)
                    mind = min (min (min (c1, c2), c3), c4)
                    if mind == c4:
                       clip = r4
                    elif mind == c3:
                       clip = r3
                    elif mind == c2:
                       clip = r2
                    else:
                       clip = r1
                    dst_plane[y][x] = clip
        return fout

    def rg7 (n, f):
        fout = f.copy ()
    
        for p in range (fout.format.num_planes):
            plane        = f.get_read_array (p)
            dst_plane    = fout.get_write_array (p)
    
            plane_height = len (plane)
            plane_width  = len (plane[0])
     
            for y in range (1, plane_height - 1):
                for x in range (1, plane_width - 1):
                    max1 = max (plane[y-1][x-1], plane[y+1][x+1])
                    min1 = min (plane[y-1][x-1], plane[y+1][x+1])
                    max2 = max (plane[y-1][x], plane[y+1][x])
                    min2 = min (plane[y-1][x], plane[y+1][x])
                    max3 = max (plane[y-1][x+1], plane[y+1][x-1])
                    min3 = min (plane[y-1][x+1], plane[y+1][x-1])
                    max4 = max (plane[y][x-1], plane[y][x+1])
                    min4 = min (plane[y][x-1], plane[y][x+1])
                    d1   = max1 - min1
                    d2   = max2 - min2
                    d3   = max3 - min3
                    d4   = max4 - min4
                    r4   = max (min (plane[y][x], max4), min4)
                    r3   = max (min (plane[y][x], max3), min3)
                    r2   = max (min (plane[y][x], max2), min2)
                    r1   = max (min (plane[y][x], max1), min1)
                    c1   = abs (plane[y][x] - r1) + d1
                    c2   = abs (plane[y][x] - r2) + d2
                    c3   = abs (plane[y][x] - r3) + d3
                    c4   = abs (plane[y][x] - r4) + d4
                    mind = min (min (min (c1, c2), c3), c4)
                    if mind == c4:
                       clip = r4
                    elif mind == c3:
                       clip = r3
                    elif mind == c2:
                       clip = r2
                    else:
                       clip = r1
                    dst_plane[y][x] = clip
        return fout

    def rg8 (n, f):
        fout = f.copy ()
    
        for p in range (fout.format.num_planes):
            plane        = f.get_read_array (p)
            dst_plane    = fout.get_write_array (p)
    
            plane_height = len (plane)
            plane_width  = len (plane[0])
     
            for y in range (1, plane_height - 1):
                for x in range (1, plane_width - 1):
                    max1 = max (plane[y-1][x-1], plane[y+1][x+1])
                    min1 = min (plane[y-1][x-1], plane[y+1][x+1])
                    max2 = max (plane[y-1][x], plane[y+1][x])
                    min2 = min (plane[y-1][x], plane[y+1][x])
                    max3 = max (plane[y-1][x+1], plane[y+1][x-1])
                    min3 = min (plane[y-1][x+1], plane[y+1][x-1])
                    max4 = max (plane[y][x-1], plane[y][x+1])
                    min4 = min (plane[y][x-1], plane[y][x+1])
                    d1   = max1 - min1
                    d2   = max2 - min2
                    d3   = max3 - min3
                    d4   = max4 - min4
                    r4   = max (min (plane[y][x], max4), min4)
                    r3   = max (min (plane[y][x], max3), min3)
                    r2   = max (min (plane[y][x], max2), min2)
                    r1   = max (min (plane[y][x], max1), min1)
                    c1   = max (min (abs (plane[y][x] - r1) + (d1 * 2), 1.0), 0.0)
                    c2   = max (min (abs (plane[y][x] - r2) + (d2 * 2), 1.0), 0.0)
                    c3   = max (min (abs (plane[y][x] - r3) + (d3 * 2), 1.0), 0.0)
                    c4   = max (min (abs (plane[y][x] - r4) + (d4 * 2), 1.0), 0.0)
                    mind = min (min (min (c1, c2), c3), c4)
                    if mind == c4:
                       clip = r4
                    elif mind == c3:
                       clip = r3
                    elif mind == c2:
                       clip = r2
                    else:
                       clip = r1
                    dst_plane[y][x] = clip
        return fout

    def rg9 (n, f):
        fout = f.copy ()
    
        for p in range (fout.format.num_planes):
            plane        = f.get_read_array (p)
            dst_plane    = fout.get_write_array (p)
    
            plane_height = len (plane)
            plane_width  = len (plane[0])
     
            for y in range (1, plane_height - 1):
                for x in range (1, plane_width - 1):
                    max1 = max (plane[y-1][x-1], plane[y+1][x+1])
                    min1 = min (plane[y-1][x-1], plane[y+1][x+1])
                    max2 = max (plane[y-1][x], plane[y+1][x])
                    min2 = min (plane[y-1][x], plane[y+1][x])
                    max3 = max (plane[y-1][x+1], plane[y+1][x-1])
                    min3 = min (plane[y-1][x+1], plane[y+1][x-1])
                    max4 = max (plane[y][x-1], plane[y][x+1])
                    min4 = min (plane[y][x-1], plane[y][x+1])
                    d1   = max1 - min1
                    d2   = max2 - min2
                    d3   = max3 - min3
                    d4   = max4 - min4
                    r4   = max (min (plane[y][x], max4), min4)
                    r3   = max (min (plane[y][x], max3), min3)
                    r2   = max (min (plane[y][x], max2), min2)
                    r1   = max (min (plane[y][x], max1), min1)
                    mind = min (min (min (d1, d2), d3), d4)
                    if mind == d4:
                       clip = r4
                    elif mind == d3:
                       clip = r3
                    elif mind == d2:
                       clip = r2
                    else:
                       clip = r1
                    dst_plane[y][x] = clip
        return fout

    def rg10 (n, f):
        fout = f.copy ()
    
        for p in range (fout.format.num_planes):
            plane        = f.get_read_array (p)
            dst_plane    = fout.get_write_array (p)
    
            plane_height = len (plane)
            plane_width  = len (plane[0])
     
            for y in range (1, plane_height - 1):
                for x in range (1, plane_width - 1):
                    d1   = abs (plane[y][x] - plane[y-1][x-1])
                    d2   = abs (plane[y][x] - plane[y-1][x])
                    d3   = abs (plane[y][x] - plane[y-1][x+1])
                    d4   = abs (plane[y][x] - plane[y][x-1])
                    d5   = abs (plane[y][x] - plane[y][x+1])
                    d6   = abs (plane[y][x] - plane[y+1][x-1])
                    d7   = abs (plane[y][x] - plane[y+1][x])
                    d8   = abs (plane[y][x] - plane[y+1][x+1])
                    mind = min (min (min (min (min (min (min (d1, d2), d3), d4), d5), d6), d7), d8)
                    if mind == d8:
                       clip = plane[y+1][x+1]
                    elif mind == d7:
                       clip = plane[y+1][x]
                    elif mind == d6:
                       clip = plane[y+1][x-1]
                    elif mind == d5:
                       clip = plane[y][x+1]
                    elif mind == d4:
                       clip = plane[y][x-1]
                    elif mind == d3:
                       clip = plane[y-1][x+1]
                    elif mind == d2:
                       clip = plane[y-1][x]
                    else:
                       clip = plane[y-1][x-1]
                    dst_plane[y][x] = clip
        return fout

    def rg11 (n, f):
        fout = f.copy ()
    
        for p in range (fout.format.num_planes):
            plane        = f.get_read_array (p)
            dst_plane    = fout.get_write_array (p)
    
            plane_height = len (plane)
            plane_width  = len (plane[0])
       
            for y in range (1, plane_height - 1):
                for x in range (1, plane_width - 1):
                    dst_plane[y][x] = (plane[y][x] * 4 + (plane[y+1][x] + plane[y-1][x] + plane[y][x+1] + plane[y][x-1]) * 2 + plane[y+1][x+1] + plane[y+1][x-1] + plane[y-1][x+1] + plane[y-1][x-1]) / 16
        return fout

    def rg13 (n, f):
        fout = f.copy ()
    
        for p in range (fout.format.num_planes):
            plane        = f.get_read_array (p)
            dst_plane    = fout.get_write_array (p)
    
            plane_height = len (plane)
            plane_width  = len (plane[0])
     
            for y in range (1, plane_height - 1):
                for x in range (1, plane_width - 1):
                    d1   = abs (plane[y+1][x-1] - plane[y-1][x+1])
                    d2   = abs (plane[y+1][x] - plane[y-1][x])
                    d3   = abs (plane[y+1][x+1] - plane[y-1][x-1])
                    mind = min (min (d1, d2), d3)
                    if mind == d3:
                       clip = (plane[y+1][x+1] + plane[y-1][x-1]) / 2
                    elif mind == d2:
                       clip = (plane[y+1][x] + plane[y-1][x]) / 2
                    else:
                       clip = (plane[y+1][x-1] + plane[y-1][x+1]) / 2
                    dst_plane[y][x] = clip
        return fout

    def rg15 (n, f):
        fout = f.copy ()
    
        for p in range (fout.format.num_planes):
            plane        = f.get_read_array (p)
            dst_plane    = fout.get_write_array (p)
    
            plane_height = len (plane)
            plane_width  = len (plane[0])
     
            for y in range (1, plane_height - 1):
                for x in range (1, plane_width - 1):
                    d1   = abs (plane[y+1][x-1] - plane[y-1][x+1])
                    d2   = abs (plane[y+1][x] - plane[y-1][x])
                    d3   = abs (plane[y+1][x+1] - plane[y-1][x-1])
                    min1 = min (plane[y+1][x-1], plane[y-1][x+1])
                    max1 = max (plane[y+1][x-1], plane[y-1][x+1])
                    min2 = min (plane[y+1][x], plane[y-1][x])
                    max2 = max (plane[y+1][x], plane[y-1][x])
                    min3 = min (plane[y+1][x+1], plane[y-1][x-1])
                    max3 = max (plane[y+1][x+1], plane[y-1][x-1])
                    mind = min (min (d1, d2), d3)
                    avrg = (plane[y+1][x-1] + plane[y+1][x] * 2 + plane[y+1][x+1] + plane[y-1][x-1] + plane[y-1][x] * 2 + plane[y-1][x+1]) / 8
                    if mind == d3:
                       clip = max (min (avrg, max3), min3)
                    elif mind == d2:
                       clip = max (min (avrg, max2), min2)
                    else:
                       clip = max (min (avrg, max1), min1)
                    dst_plane[y][x] = clip
        return fout

    def rg17 (n, f):
        fout = f.copy ()
    
        for p in range (fout.format.num_planes):
            plane        = f.get_read_array (p)
            dst_plane    = fout.get_write_array (p)
    
            plane_height = len (plane)
            plane_width  = len (plane[0])
     
            for y in range (1, plane_height - 1):
                for x in range (1, plane_width - 1):
                    max1 = max (plane[y-1][x-1], plane[y+1][x+1])
                    min1 = min (plane[y-1][x-1], plane[y+1][x+1])
                    max2 = max (plane[y-1][x], plane[y+1][x])
                    min2 = min (plane[y-1][x], plane[y+1][x])
                    max3 = max (plane[y-1][x+1], plane[y+1][x-1])
                    min3 = min (plane[y-1][x+1], plane[y+1][x-1])
                    max4 = max (plane[y][x-1], plane[y][x+1])
                    min4 = min (plane[y][x-1], plane[y][x+1])
                    low  = max (max (max (min1, min2), min3), min4)
                    up   = min (min (min (max1, max2), max3), max4)
                    maxn = max (low, up)
                    minn = min (low, up)
                    dst_plane[y][x] = max (min (plane[y][x], maxn), minn)
        return fout

    def rg18 (n, f):
        fout = f.copy ()
    
        for p in range (fout.format.num_planes):
            plane        = f.get_read_array (p)
            dst_plane    = fout.get_write_array (p)
    
            plane_height = len (plane)
            plane_width  = len (plane[0])
     
            for y in range (1, plane_height - 1):
                for x in range (1, plane_width - 1):
                    d1   = max (abs (plane[y][x] - plane[y-1][x-1]), abs (plane[y][x] - plane[y+1][x+1]))
                    d2   = max (abs (plane[y][x] - plane[y-1][x]), abs (plane[y][x] - plane[y+1][x]))
                    d3   = max (abs (plane[y][x] - plane[y-1][x+1]), abs (plane[y][x] - plane[y+1][x-1]))
                    d4   = max (abs (plane[y][x] - plane[y][x-1]), abs (plane[y][x] - plane[y][x+1]))
                    max1 = max (plane[y-1][x-1], plane[y+1][x+1])
                    min1 = min (plane[y-1][x-1], plane[y+1][x+1])
                    max2 = max (plane[y-1][x], plane[y+1][x])
                    min2 = min (plane[y-1][x], plane[y+1][x])
                    max3 = max (plane[y-1][x+1], plane[y+1][x-1])
                    min3 = min (plane[y-1][x+1], plane[y+1][x-1])
                    max4 = max (plane[y][x-1], plane[y][x+1])
                    min4 = min (plane[y][x-1], plane[y][x+1])
                    r4   = max (min (plane[y][x], max4), min4)
                    r3   = max (min (plane[y][x], max3), min3)
                    r2   = max (min (plane[y][x], max2), min2)
                    r1   = max (min (plane[y][x], max1), min1)
                    mind = min (min (min (d1, d2), d3), d4)
                    if mind == d4:
                       clip = r4
                    elif mind == d3:
                       clip = r3
                    elif mind == d2:
                       clip = r2
                    else:
                       clip = r1
                    dst_plane[y][x] = clip
        return fout

    def rg19 (n, f):
        fout = f.copy ()
    
        for p in range (fout.format.num_planes):
            plane        = f.get_read_array (p)
            dst_plane    = fout.get_write_array (p)
    
            plane_height = len (plane)
            plane_width  = len (plane[0])
       
            for y in range (1, plane_height - 1):
                for x in range (1, plane_width - 1):
                    dst_plane[y][x] = (plane[y+1][x] + plane[y-1][x] + plane[y][x+1] + plane[y][x-1] + plane[y+1][x+1] + plane[y+1][x-1] + plane[y-1][x+1] + plane[y-1][x-1]) / 8
        return fout

    def rg20 (n, f):
        fout = f.copy ()
    
        for p in range (fout.format.num_planes):
            plane        = f.get_read_array (p)
            dst_plane    = fout.get_write_array (p)
    
            plane_height = len (plane)
            plane_width  = len (plane[0])
       
            for y in range (1, plane_height - 1):
                for x in range (1, plane_width - 1):
                    dst_plane[y][x] = (plane[y][x] + plane[y+1][x] + plane[y-1][x] + plane[y][x+1] + plane[y][x-1] + plane[y+1][x+1] + plane[y+1][x-1] + plane[y-1][x+1] + plane[y-1][x-1]) / 9
        return fout

    def rg21 (n, f):
        fout = f.copy ()
    
        for p in range (fout.format.num_planes):
            plane        = f.get_read_array (p)
            dst_plane    = fout.get_write_array (p)
    
            plane_height = len (plane)
            plane_width  = len (plane[0])
     
            for y in range (1, plane_height - 1):
                for x in range (1, plane_width - 1):
                    avg1 = (plane[y-1][x-1] + plane[y+1][x+1]) / 2
                    avg2 = (plane[y-1][x] + plane[y+1][x]) / 2
                    avg3 = (plane[y-1][x+1] + plane[y+1][x-1]) / 2
                    avg4 = (plane[y][x-1] + plane[y][x+1]) / 2
                    maxn = max (max (max (avg1, avg2), avg3), avg4)
                    minn = min (min (min (avg1, avg2), avg3), avg4)
                    dst_plane[y][x] = max (min (plane[y][x], maxn), minn)
        return fout

    def rg23 (n, f):
        fout = f.copy ()
    
        for p in range (fout.format.num_planes):
            plane        = f.get_read_array (p)
            dst_plane    = fout.get_write_array (p)
    
            plane_height = len (plane)
            plane_width  = len (plane[0])
     
            for y in range (1, plane_height - 1):
                for x in range (1, plane_width - 1):
                    max1 = max (plane[y-1][x-1], plane[y+1][x+1])
                    min1 = min (plane[y-1][x-1], plane[y+1][x+1])
                    max2 = max (plane[y-1][x], plane[y+1][x])
                    min2 = min (plane[y-1][x], plane[y+1][x])
                    max3 = max (plane[y-1][x+1], plane[y+1][x-1])
                    min3 = min (plane[y-1][x+1], plane[y+1][x-1])
                    max4 = max (plane[y][x-1], plane[y][x+1])
                    min4 = min (plane[y][x-1], plane[y][x+1])
                    d1   = max1 - min1
                    d2   = max2 - min2
                    d3   = max3 - min3
                    d4   = max4 - min4
                    u1   = min (plane[y][x] - max1, d1)
                    u2   = min (plane[y][x] - max2, d2)
                    u3   = min (plane[y][x] - max3, d3)
                    u4   = min (plane[y][x] - max4, d4)
                    u    = max (max (max (max (u1, u2), u3), u4), 0.0)
                    l1   = min (min1 - plane[y][x], d1)
                    l2   = min (min2 - plane[y][x], d2)
                    l3   = min (min3 - plane[y][x], d3)
                    l4   = min (min4 - plane[y][x], d4)
                    l    = max (max (max (max (l1, l2), l3), l4), 0.0)
                    dst_plane[y][x] = plane[y][x] - u + l
        return fout

    def rg24 (n, f):
        fout = f.copy ()
    
        for p in range (fout.format.num_planes):
            plane        = f.get_read_array (p)
            dst_plane    = fout.get_write_array (p)
    
            plane_height = len (plane)
            plane_width  = len (plane[0])
     
            for y in range (1, plane_height - 1):
                for x in range (1, plane_width - 1):
                    max1 = max (plane[y-1][x-1], plane[y+1][x+1])
                    min1 = min (plane[y-1][x-1], plane[y+1][x+1])
                    max2 = max (plane[y-1][x], plane[y+1][x])
                    min2 = min (plane[y-1][x], plane[y+1][x])
                    max3 = max (plane[y-1][x+1], plane[y+1][x-1])
                    min3 = min (plane[y-1][x+1], plane[y+1][x-1])
                    max4 = max (plane[y][x-1], plane[y][x+1])
                    min4 = min (plane[y][x-1], plane[y][x+1])
                    d1   = max1 - min1
                    d2   = max2 - min2
                    d3   = max3 - min3
                    d4   = max4 - min4
                    t1   = plane[y][x] - max1
                    t2   = plane[y][x] - max2
                    t3   = plane[y][x] - max3
                    t4   = plane[y][x] - max4
                    u1   = min (t1, d1 - t1)
                    u2   = min (t2, d2 - t2)
                    u3   = min (t3, d3 - t3)
                    u4   = min (t4, d4 - t4)
                    u    = max (max (max (max (u1, u2), u3), u4), 0.0)
                    t1   = min1 - plane[y][x]
                    t2   = min2 - plane[y][x]
                    t3   = min3 - plane[y][x]
                    t4   = min4 - plane[y][x]
                    l1   = min (t1, d1 - t1)
                    l2   = min (t2, d2 - t2)
                    l3   = min (t3, d3 - t3)
                    l4   = min (t4, d4 - t4)
                    l    = max (max (max (max (l1, l2), l3), l4), 0.0)
                    dst_plane[y][x] = plane[y][x] - u + l
        return fout
    if mode == 1:
       mode = rg1
    elif mode == 2:
       mode = rg2
    elif mode == 3:
       mode = rg3
    elif mode == 4:
       mode = rg4
    elif mode == 5:
       mode = rg5
    elif mode == 6:
       mode = rg6
    elif mode == 7:
       mode = rg7
    elif mode == 8:
       mode = rg8
    elif mode == 9:
       mode = rg9
    elif mode == 10:
       mode = rg10
    elif mode == 11:
       mode = rg11
    elif mode == 12:
       mode = rg11
    elif mode == 13:
       mode = rg13
    elif mode == 14:
       mode = rg13
    elif mode == 15:
       mode = rg15
    elif mode == 16:
       mode = rg15
    elif mode == 17:
       mode = rg17
    elif mode == 18:
       mode = rg18
    elif mode == 19:
       mode = rg19
    elif mode == 20:
       mode = rg20
    elif mode == 21:
       mode = rg21
    elif mode == 22:
       mode = rg21
    elif mode == 23:
       mode = rg23
    elif mode == 24:
       mode = rg24
    else:
       mode = rg0
    clip = padding (src, 1, 1, 1, 1)
    clip = core.std.ModifyFrame (clip=clip, clips=clip, selector=mode)
    clip = core.std.CropRel (clip, 1, 1, 1, 1)
    return clip
