# Compute matrix ratio using serpentine binning
# cmdoret, 20190919

import sys
import numpy as np
import hicstuff.hicstuff as hcs
import hicstuff.view as hcv
import hicstuff.io as hio
import serpentine as sp
import matplotlib.pyplot as plt

A = sys.argv[1]
B = sys.argv[2]
out = sys.argv[3]

# Load matrix data from cool
A, fa, ca = hio.load_cool(A)
B, fb, cb = hio.load_cool(B)

# Normalise input matrices
A = hcv.sparse_to_dense(A, remove_diag=False)
B = hcv.sparse_to_dense(B, remove_diag=False)

# Filter out low detectability bins
flt = sp.outstanding_filter(A) + sp.outstanding_filter(B)
flt = flt == False
A = sp.fltmatr(A, flt)
B = sp.fltmatr(B, flt)

# Find the binning threshold
mdtrend, mdthreshold = sp.MDbefore(A, B, show=False, s=10)

threshold = mdthreshold
minthreshold = threshold / 5.0

sA, sB, sK = sp.serpentin_binning(
    A,
    B,
    threshold=threshold,
    minthreshold=minthreshold,
    triangular=True,
    parallel=2,
)

sp.serpentine._plot(sA, sB, sK - mdtrend, triangular=False, limit=3)
plt.savefig(out)
