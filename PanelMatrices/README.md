# PanelMatrices

`PanelMatrices` is work in progress. It will provide `Matrix` types that use panel-major storage, and utilize hand-written assembly kernels from BLASFEO for performance.

# Plan

* A `PanelMatrix` will contain its size and panel type as type parameters. (All of the parameters needed by BLASFEO will be stored in the type - not in a struct.)

* Generic Julia code will be written to handle the panel operations.

* For certain (`eltype`,`size`) combinations, specialized kernels from BLASFEO will override the generic code.
