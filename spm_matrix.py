import numpy as np
import scipy
import matcompat

# if available import pylab (from matlibplot)
try:
    import matplotlib.pylab as plt
except ImportError:
    pass


def spm_matrix(P, order):
    # Local Variables: A, R1, R2, R3, q, P, S, R, T, Z, order
    # Function calls: cos, eye, eval, error, nargin, length, spm_matrix, sprintf, numel, isequal, isnumeric, sin, size
    # % Return an affine transformation matrix
    # % FORMAT [A] = spm_matrix(P [,order])
    # % P(1)  - x translation
    # % P(2)  - y translation
    # % P(3)  - z translation
    # % P(4)  - x rotation about - {pitch} (radians)
    # % P(5)  - y rotation about - {roll}  (radians)
    # % P(6)  - z rotation about - {yaw}   (radians)
    # % P(7)  - x scaling
    # % P(8)  - y scaling
    # % P(9)  - z scaling
    # % P(10) - x affine
    # % P(11) - y affine
    # % P(12) - z affine
    # %
    # % order - application order of transformations [Default: 'T*R*Z*S']
    # %
    # % A     - affine transformation matrix
    # %__________________________________________________________________________
    # %
    # % spm_matrix returns a matrix defining an orthogonal linear (translation,
    # % rotation, scaling or affine) transformation given a vector of
    # % parameters (P).  By default, the transformations are applied in the
    # % following order (i.e., the opposite to which they are specified):
    # %
    # % 1) shear
    # % 2) scale (zoom)
    # % 3) rotation - yaw, roll & pitch
    # % 4) translation
    # %
    # % This order can be changed by calling spm_matrix with a string as a
    # % second argument. This string may contain any valid MATLAB expression
    # % that returns a 4x4 matrix after evaluation. The special characters 'S',
    # % 'Z', 'R', 'T' can be used to reference the transformations 1)-4)
    # % above. The default order is 'T*R*Z*S', as described above.
    # %
    # % SPM uses a PRE-multiplication format i.e. Y = A*X where X and Y are 4 x n
    # % matrices of n coordinates.
    # %__________________________________________________________________________
    # %
    # % See also: spm_imatrix.m
    # %__________________________________________________________________________
    # % Copyright (C) 1994-2011 Wellcome Trust Centre for Neuroimaging
    # % Karl Friston
    # % $Id: spm_matrix.m 4414 2011-08-01 17:51:40Z guillaume $
    # %-Special case: translation only
    # %--------------------------------------------------------------------------
    if len(P) == 3:
        A = np.eye(4)
        A[0:3, 3] = P.flatten(1)
        return []

    # %-Pad P with 'null' parameters
    # %--------------------------------------------------------------------------
    q = np.array(np.hstack((0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0)))
    P = np.array(np.hstack((P, q[int(len(P) + 1) - 1:12])))
    # %-Translation / Rotation / Scale / Shear
    # %--------------------------------------------------------------------------
    T = np.array(np.vstack((np.hstack((1, 0, 0, P[0])), np.hstack((0, 1, 0, P[1])), np.hstack((0, 0, 1, P[2])),
                            np.hstack((0, 0, 0, 1)))))
    R1 = np.array(np.vstack((np.hstack((1, 0, 0, 0)), np.hstack((0, np.cos(P[3]), np.sin(P[3]), 0)),
                             np.hstack((0, -np.sin(P[3]), np.cos(P[3]), 0)), np.hstack((0, 0, 0, 1)))))
    R2 = np.array(np.vstack((np.hstack((np.cos(P[4]), 0, np.sin(P[4]), 0)), np.hstack((0, 1, 0, 0)),
                             np.hstack((-np.sin(P[4]), 0, np.cos(P[4]), 0)), np.hstack((0, 0, 0, 1)))))
    R3 = np.array(np.vstack((np.hstack((np.cos(P[5]), np.sin(P[5]), 0, 0)),
                             np.hstack((-np.sin(P[5]), np.cos(P[5]), 0, 0)), np.hstack((0, 0, 1, 0)),
                             np.hstack((0, 0, 0, 1)))))
    R = np.dot(np.dot(R1, R2), R3)
    Z = np.array(np.vstack((np.hstack((P[6], 0, 0, 0)), np.hstack((0, P[7], 0, 0)), np.hstack((0, 0, P[8], 0)),
                            np.hstack((0, 0, 0, 1)))))
    S = np.array(np.vstack((np.hstack((1, P[9], P[10], 0)), np.hstack((0, 1, P[11], 0)),
                            np.hstack((0, 0, 1, 0)), np.hstack((0, 0, 0, 1)))))
    # %-Affine transformation matrix
    # %--------------------------------------------------------------------------

    A = np.dot(np.dot(np.dot(T, R), Z), S)

    return [A]
