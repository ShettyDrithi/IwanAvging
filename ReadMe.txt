Averaging algorithm for integrating SDOF(/modal) system containing an Iwan element forced in the microslip regime

-->iwanavg_example.m is an example script that calls the main function iwanavg.m

[resp,simtime,varargout] = iwanavg(Forcing,time,paramp,m,Kinf,zt_lin,varargin)

--> required input for iwanavg.m:
	Force, either as a function of time or a vector. if no forcing present, enter a dummy zero-amplitude force
	(need to set up for no forcing)
	time, either simulation time or a vector that has same spacing as the force vector. when only simulation time is entered,
	an equally spaced time vector with dt=1/(fn0*200) is created. spacing chosen arbitrarily based on what worked best for
	numerical schemes.
	physical Iwan parameters [Fs,Kt,chi,beta]
	mass
	linear, macroslip stiffness (Kinf) and material damping (ztlin)
		optional input: enter as name-value pairs (NOT case sensitive)
		TBand: The time band to fit piecewise-hilbert transform. Hilbert not performed if Tband not provided
		ICs: Initial conditions
--> output of iwanavg.m:
	resp: contains four rows - displacement, velocity, accn and time vectors
	simtime: total time taken to complete simulation
		optional output:
		wn_fit,zt_fit,yfit - frequency and damping vs velocity amplitude obtained from Hilbert
		Ya,Yv - FFTs of accn and velocity vs frequency ws in rad/s
		wn, zeta - frequency and damping obtained directly from averaging (without Hilbert) with corresponding 
			   time vector tfree. you can also use velocity vector vfree but that will require some post-processing
		           to find the peaks/amplitude
