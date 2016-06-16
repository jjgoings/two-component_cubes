This is a python class to read and write out Gaussian cube files, 
with a particular emphasisis on extracting complex and two-component molecular 
orbitals.

Example:
Say you have a complex 2c MO cube called `twoc.cube`.
If you have `cube.py` in the same directory, you could write a script like
```python
from cube import Cube
twoc = Cube('twoc.cube')
twoc.write_out('twoc_ra.cube',data='RA')
twoc.write_out('twoc_ia.cube',data='IA')
twoc.write_out('twoc_rb.cube',data='RB')
twoc.write_out('twoc_ib.cube',data='IB')
```

This would create an object called `twoc`, which you could then dump out new
cube files that are readable by GaussView corresponding to the real alpha
(`data='RA'`), real beta (`data='RB'`), imaginary alpha (`data='IA'`), and
imaginary beta (`data='IB'`). The default for `write_out` is real alpha.  

The code should be smart enough to distinguish between MO and density cube
files.

Usually there is no problem using GaussView for RHF and UHF anyway, so the code
isn't yet generalized to all cases, but I'll work that out when I have time. 

 
