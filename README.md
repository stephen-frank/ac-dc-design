# AC-DC Design Source Code

This electronic appendix includes the GAMS source code for the formulation and
NGBD solution algorithm for the mixed AC-DC electricity distribution system
design problem described in the article:

> S. M. Frank and S. Rebennack, "Optimal Design of Mixed AC-DC Distribution
> Systems for Commercial Buildings: A Nonconvex Generalized Benders
> Decomposition Approach," *European Journal of Operational Research*, 2014,
> accepted for publication.

Much of the source code is extended from the following dissertation:
    
> S. M. Frank, "Optimal Design of Mixed AC-DC Distribution Systems for
> Commercial Buildings," Dissertation, Colorado School of Mines, Golden, CO,
> 2013. Available: http://disexpress.umi.com/dxweb, UMI #3558161.
    
In addition to the descriptions provided in the article and dissertation, the
code in this appendix is self-documented via extensive comments at the head of
and within each file. The `Test`, `Batch`, and `Case_Studies` directories also
include local README files which describe the directory contents and provide
some usage guidance. If you're looking for a place to start, try the README in
the `Test` directory.

**Please note that the underlying software products I used to create this code
(MATLAB and GAMS) have changed since I completed my dissertation. Therefore,
this code no longer runs "out-of-the-box". Please consider this repository a
reference and code archive rather than working software.** See also the
**Notes** section below.
    
This is version 1.1 of the electronic appendix. The latest version of this
appendix may be found (after article publication) at Stephen Frank's personal
website, http://www.stevefrank.info or on Github at
https://github.com/stephen-frank/ac-dc-design/. This code is released under
the Gnu Public License (GPL) version 3. A copy of this license is available
in the root directory of this appendix (the same directory as this README file).

## Notes
1. This MATLAB code in this repository was created for MATLAB 2013. Users have reported
that subsequent versions of MATLAB (2014 and later) are not compatible with
the object-oriented approach I implemented. Therefore, the code will not run
as-is with later MATLAB versions.
