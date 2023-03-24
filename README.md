# 001_SteadyHeatdiffusion_with_SourceTerm
Usual starting point: 1D steady heat diffusion equation with a source term, with prescribed temperature on either ends of the domain.
1. The program assumes a 1D domain with length L. 
2. Prescribed temperature boundary conditions on either side of the domain.
3. Ghost points used to specify boundary condition on the RHS matrix: for example; (1/2) * (T-1 + T0) = TL
4. Source term is included in RHS matrix. 
5. At the moment: A source term can be added using direct array assignment --> for example S[4] = A 

