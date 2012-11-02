import sys
from math import sin, cos

identity = [1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0]

def mult(A, B):
    return [
        A[0]*B[0]+A[4]*B[1]+A[8]*B[2]+A[12]*B[3],
        A[1]*B[0]+A[5]*B[1]+A[9]*B[2]+A[13]*B[3],
        A[2]*B[0]+A[6]*B[1]+A[10]*B[2]+A[14]*B[3],
        A[3]*B[0]+A[7]*B[1]+A[11]*B[2]+A[15]*B[3],
        A[0]*B[4]+A[4]*B[5]+A[8]*B[6]+A[12]*B[7],
        A[1]*B[4]+A[5]*B[5]+A[9]*B[6]+A[13]*B[7],
        A[2]*B[4]+A[6]*B[5]+A[10]*B[6]+A[14]*B[7],
        A[3]*B[4]+A[7]*B[5]+A[11]*B[6]+A[15]*B[7],
        A[0]*B[8]+A[4]*B[9]+A[8]*B[10]+A[12]*B[11],
        A[1]*B[8]+A[5]*B[9]+A[9]*B[10]+A[13]*B[11],
        A[2]*B[8]+A[6]*B[9]+A[10]*B[10]+A[14]*B[11],
        A[3]*B[8]+A[7]*B[9]+A[11]*B[10]+A[15]*B[11],
        A[0]*B[12]+A[4]*B[13]+A[8]*B[14]+A[12]*B[15],
        A[1]*B[12]+A[5]*B[13]+A[9]*B[14]+A[13]*B[15],
        A[2]*B[12]+A[6]*B[13]+A[10]*B[14]+A[14]*B[15],
        A[3]*B[12]+A[7]*B[13]+A[11]*B[14]+A[15]*B[15]
            ]

def get_rot_x(theta):
    return [1.0, 0.0, 0.0, 0.0, 0.0, cos(theta), sin(theta), 0.0, 0.0, -sin(theta), cos(theta), 0.0, 0.0, 0.0, 0.0, 1.0]

def get_rot_y(theta):
    return [cos(theta), 0.0, -sin(theta), 0.0, 0.0, 1.0, 0.0, 0.0, sin(theta), 0.0, cos(theta), 0.0, 0.0, 0.0, 0.0, 1.0]

def get_rot_z(theta):
    return [cos(theta), sin(theta), 0.0, 0.0, -sin(theta), cos(theta), 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0]

def rotate_x(A, theta):
    return mult(get_rot_x(theta),A)

def rotate_y(A, theta):
    return mult(get_rot_y(theta),A)

def rotate_z(A, theta):
    return mult(get_rot_z(theta),A)

def transl(A, x, y, z):
    return [A[0], A[1], A[2], A[3], A[4], A[5], A[6], A[7], A[8], A[9], A[10], A[11], A[12]+x, A[13]+y, A[14]+z, A[15]]

if __name__ == "__main__":
    # get the final transformation
    with open(sys.argv[1], "r") as f:
        for line in f:
            mat = line.split()
            if len(mat) != 17:
                print "invalid line: %s"%line
                exit(1)
            mat = [float(i) for i in mat[:16]] # discard color

    x, y, z = mat[12:15]
    res = []

    rotdelta = 0.001
    rot = [0]
    for f in range(100):
        rotdelta*=1.07
        rot.append(rot[-1]+rotdelta)
    for r in reversed(rot):
        m = transl(mat, -x, -y, -z)
        m = rotate_y(m, r)
        m = transl(m, x, y, z)
        res.append(m)

    with open(sys.argv[1], "w") as f:
        for m in res:
            f.write("%s\n"%" ".join([str(i) for i in m+[1]]))

if __name__ == "foo":
    with open(sys.argv[1]) as f:
        for line in f:
            mat = line.split()
            if len(mat) != 17:
                print "invalid line: %s"%line
                exit(1)
            mat = [float(i) for i in mat[:16]] # discard color
            x, y, z = mat[12:15]
            mat = transl(mat, -x, -y, -z)
            mat = rotate_y(mat, rot)
            mat = transl(mat, x, y, z)
            print " ".join([str(i) for i in mat+[1]])
            rot += rotdelta
