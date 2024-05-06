
using LaMEM, GeophysicalModelGenerator, Plots
include("subindshapes.jl")  # Include the file containing shape definitions and functions
include("custom_gmg.jl")  # Include the file containing custom functions for GeophysicalModelGenerator

modelHeight = 2880.0km
plateHeight = 120.0km
T_mantle = 1350  # in Celcius

refDensity = 3200.0kg / m^3
deltaRhoMax = 80.0kg / m^3
gravity = 9.81m / s^2
refViscosity = 5e20Pas
bodyForce = refDensity * gravity
K_tau = bodyForce * modelHeight

CharUnits = GEO_units(length=modelHeight, stress=K_tau, viscosity=refViscosity)

modelHeight = modelHeight.val

resFactor = 1.0#rf

model = Model(
    Grid(
        coord_x=[0, 5040, 6192, 11520],
        nel_x=round.(Int, [365, 480, 435] / resFactor),
        # nel_x=round.(Int, [365, 480*2, 435] / resFactor),

        bias_x=[0.125, 1.0, 4],
        y=[-2.5, 2.5],
        nel_y=1,
        coord_z=[-2880, -384, 0, 192],
        # coord_z=[-2880, -384, 0],
        nel_z=round.(Int, [320, 256, 64] / resFactor),
        bias_z=[0.125, 1.0, 4]
    ),
    BoundaryConditions(open_top_bound=1), # (left right front back bottom top)
    Scaling(CharUnits)
)

refViscosity = 5e20
opViscFactor = 1.0
modelMaterials = [
    Phase(ID=0, Name="air", rho=50.0, eta=1e-2 * refViscosity),
    Phase(ID=1, Name="UpperMantle", eta=refViscosity, eta0=refViscosity, n=3.5, e0=2.5e-15, rho=3200.0, eta_st=0.1 * refViscosity, eta_vp=0.1 * refViscosity),
    Phase(ID=2, Name="UppperCrustIndoAustralianPlate", eta=1e3 * refViscosity, rho=3280.0, ch=6e6, eta_vp=0.1 * refViscosity, eta_st=0.1 * refViscosity),
    Phase(ID=3, Name="LowerCrustIndoAustralianPlate", eta=1e3 * refViscosity, rho=3280.0, ch=50e6, eta_vp=refViscosity, eta_st=refViscosity),
    Phase(ID=4, Name="LithosphericMantleCrustIndoAustralianPlate", eta=1e3 * refViscosity, rho=3280.0, ch=350e6, eta_vp=refViscosity, eta_st=refViscosity),
    Phase(ID=5, Name="LithosphericMantleIndoAustralianPlate", eta=5e1 * refViscosity, rho=3280.0, ch=30e6, eta_vp=refViscosity, eta_st=refViscosity), Phase(ID=6, Name="UpperCrustIndianIndentor", eta=1e3 * refViscosity, rho=2800.0, ch=06e6, eta_vp=0.1 * refViscosity, eta_st=0.1 * refViscosity),
    Phase(ID=7, Name="UpperCrustIndianIndentor", eta=1e3 * refViscosity, rho=2800.0, ch=100e6, eta_vp=refViscosity, eta_st=refViscosity), Phase(ID=8, Name="WeakULCrustIndianIndentor", eta=1e3 * refViscosity, rho=3000.0, ch=06e6, eta_vp=0.1 * refViscosity, eta_st=0.1 * refViscosity),
    Phase(ID=9, Name="LowerCrustIndianIndentor", eta=1e3 * refViscosity, rho=3100.0, ch=200e6, eta_vp=refViscosity, eta_st=refViscosity),
    Phase(ID=10, Name="WeakLLCrustIndianIndentor", eta=1e3 * refViscosity, rho=3100.0, ch=06e6, eta_vp=0.1 * refViscosity, eta_st=0.1 * refViscosity), Phase(ID=11, Name="LithosphericMantleIndianIndentor", eta=1e3 * refViscosity, rho=3100.0, ch=350e6, eta_vp=refViscosity, eta_st=refViscosity),
    Phase(ID=12, Name="LithosphericMantleIndianIndentor", eta=5e1 * refViscosity, rho=3200.0, ch=30e6, eta_vp=refViscosity, eta_st=refViscosity), Phase(ID=13, Name="CrustEurasianPlateForeArc", eta=1e3 * refViscosity * opViscFactor, rho=3200.0),
    Phase(ID=14, Name="LithosphericMantleEurasianPlateForeArc", eta=5e2 * refViscosity * opViscFactor, rho=3200.0),
    Phase(ID=15, Name="CrustEurasianPlateBackArc", eta=2.5e2 * refViscosity * opViscFactor, rho=3200.0),
    Phase(ID=16, Name="LithosphericMantleEurasianPlateBackArc", eta=1.25e2 * refViscosity * opViscFactor, rho=3200.0),
    Phase(ID=17, Name="CrustEurasianPlate", eta=2e3 * refViscosity * opViscFactor, rho=3200.0),
    Phase(ID=18, Name="LithosphericMantleEurasianPlate", eta=1e3 * refViscosity * opViscFactor, rho=3200.0),
    # Phase(ID=13, Name="CrustEurasianPlateForeArc", eta=1e3 * refViscosity * opViscFactor, rho=2800.0),
    # Phase(ID=14, Name="LithosphericMantleEurasianPlateForeArc", eta=5e2 * refViscosity * opViscFactor, rho=3220.0),
    # Phase(ID=15, Name="CrustEurasianPlateBackArc", eta=2.5e2 * refViscosity * opViscFactor, rho=2800.0),
    # Phase(ID=16, Name="LithosphericMantleEurasianPlateBackArc", eta=1.25e2 * refViscosity * opViscFactor, rho=3220.0),
    # Phase(ID=17, Name="CrustEurasianPlate", eta=2e3 * refViscosity * opViscFactor, rho=2800.0),
    # Phase(ID=18, Name="LithosphericMantleEurasianPlate", eta=1e3 * refViscosity * opViscFactor, rho=3220.0),
    Phase(ID=19, Name="LowerMantle", eta=1e2 * refViscosity, rho=3200.0),]


slabshapes = Slab2D(
    top_x=2082.0, top_y=0.0, length=3672.1, taper=30.0, dip=29.0, depth=120.0,
    thickness_array=[5, 25, 30.0, 30.0])

indentorshapes = Indentor2D(
    top_x=864.0, top_y=0.0, length=2448.0, taper=18.0, taper2=12.0,
    thickness_array=[5.0, 7.5, 7.5, 12.5, 7.5, 30.0, 30.0,])

overRidingShapesForeArc = OverRidingPlate2D(
    top_x=11520.0, top_y=0.0, length=-5760.1, taper=90.0, dip=29.0,
    thickness_array=[30.0, 30.0])

overRidingShapesBackArc = OverRidingPlate2D(
    top_x=11520.0, top_y=0.0, length=-5586.0, taper=90.0, dip=90.0,
    thickness_array=[30.0, 30.0])

overRidingShapes = OverRidingPlate2D(
    top_x=11520.0, top_y=0.0, length=-5040.0, taper=90.0, dip=17.0,
    thickness_array=[40.0, 80.0])

modelshapePolygons = vcat(
    slabshapes.polygons,
    indentorshapes.polygons,
    overRidingShapesForeArc.polygons,
    overRidingShapesBackArc.polygons,
    overRidingShapes.polygons)

add_phase!(model, modelMaterials[:]...)


# Set all phases in the Grid to 1
model.Grid.Phases .= 1

# Set phases in the Grid where Z is less than -660 to 19
model.Grid.Phases[model.Grid.Grid.Z .< -660] .= 19

# Set phases in the Grid where Z is greater than 0 to 0
model.Grid.Phases[model.Grid.Grid.Z .> 0] .= 0

# Iterate over each polygon in modelshapePolygons and add it to the model
for (index, polygon) in enumerate(modelshapePolygons)
    AddPoly!(model, vertices=polygon, phase=ConstantPhase(index + 1))
end

"""
    UMLM = PhaseTransition(ID, Type, Parameter_transition, PhaseBelow, PhaseAbove, PhaseDirection, ConstantValue)

Constructs a `PhaseTransition` object with the specified parameters.

## Arguments
- `ID`: An integer representing the ID of the phase transition.
- `Type`: A string representing the type of the phase transition.
- `Parameter_transition`: A string representing the parameter transition.
- `PhaseBelow`: An array of integers representing the phases below the transition.
- `PhaseAbove`: An array of integers representing the phases above the transition.
- `PhaseDirection`: A string representing the direction of the phase transition.
- `ConstantValue`: A float representing the constant value of the phase transition.

## Returns
- `UMLM`: A `PhaseTransition` object.

"""
UMLM = PhaseTransition(
    ID=0, Type="Constant", Parameter_transition="Depth",
    PhaseBelow=[19], PhaseAbove=[1], PhaseDirection="BothWays",
    ConstantValue=-660)

# model.Materials.PhaseTransitions
add_phasetransition!(model, UMLM)

model.FreeSurface = FreeSurface(
    surf_use=1,                # free surface activation flag
    surf_corr_phase=1,                # air phase ratio correction flag (due to surface position)
    surf_level=0.0,              # initial level
    surf_air_phase=0,                # phase ID of sticky air layer
    surf_max_angle=45.0              # maximum angle with horizon (smoothed if larger))
)

model.SolutionParams = SolutionParams(
    FSSA=1,  # free surface stabilization alpha
    eta_min=1e-2 * refViscosity,
    eta_ref=1e3 * refViscosity,
    eta_max=2e3 * refViscosity,
    min_cohes=1e6,
    init_guess=1,  # initial guess flag
    act_p_shift=1,
    # DII=1e-18
)

model.Time = Time(
    time_end=150.0,
    # dt=0.001,
    dt_min=0.00000001,
    dt_max=0.1,
    nstep_max=4000,
    nstep_out=10,
    # nstep_rdb=200,
)

sloverDirect = Solver(
    SolverType="direct",
    DirectSolver="mumps",
    DirectPenalty=1e7,
    PETSc_options=[
        "-snes_ksp_ew",                     # EW scheme for the KSP
        "-snes_npicard 2",                  # Number of Picard iterations
        "-snes_monitor",                    # Monitor the nonlinear solver
        "-snes_rtol 1e-5",                  # Relative tolerance for convergence
        "-snes_atol 1e-4",                  # Absolute tolerance for convergence
        "-snes_stol 1e-16",                 # Tolerance for the step length
        "-snes_max_it 100",                  # Maximum number of nonlinear iterations
        "-snes_max_funcs 500000",           # Maximum number of function evaluations
        "-snes_max_linear_solve_fail 10000",# Maximum number of linear solve failures
        "-snes_PicardSwitchToNewton_rtol 5e-2",   # Relative tolerance for switching from Picard to Newton
        "-snes_NewtonSwitchToPicard_it 30",       # Number of Newton iterations before switching to Picard
        "-snes_NewtonSwitchToPicard_rtol 1.1",    # Relative tolerance for switching from Newton to Picard
        "-snes_linesearch_monitor",         # Monitor the line search https://petsc.org/release/manualpages/SNES/SNESNEWTONLS/
        "-snes_linesearch_type cp",         # Type of line search, cp and l2 are the most common Critical point line search. This line search assumes that there exists some artificial G(x) for which the SNESFunction F(x) = grad G(x). Therefore, this line search seeks to find roots of dot(F, Y) via a secant method.
        # "-snes_linesearch_max_it <max_it> ".
        "-snes_linesearch_maxstep 1.0",     # Maximum step length for the line search
        # "-snes_monitor_lg_residualnorm",
        "-js_ksp_type fgmres",              # Type of Krylov solver for Jacobian solve
        "-js_ksp_max_it 100",                # Maximum number of iterations for Jacobian solve
        "-js_ksp_converged_reason",         # Print the reason for convergence or divergence of the Jacobian solve
        "-js_ksp_monitor",                  # Monitor the Jacobian solve
        "-js_ksp_atol 1e-5",                # Absolute tolerance for convergence of the Jacobian solve
        "-js_ksp_rtol 1e-5",                # Relative tolerance for convergence of the Jacobian solve
        "-pcmat_type mono",                 # Type of preconditioner for the matrix
        "-pcmat_pgamma 1e3",                # Parameter for the preconditioner (similar to viscosity ratio)
        "-jp_type user",                    # Type of Jacobian preconditioner
        "-jp_pc_type lu",                   # Type of preconditioner for the Jacobian
        "-jp_pc_factor_mat_solver_type mumps",    # Solver type for the Jacobian preconditioner
        "-da_refine_y 1"                    # fancy option the grid in the y direction, 1 make it 2D they say
    ]
)

solverMultigrid = Solver(
    SolverType="multigrid",
    MGLevels=4,
    MGCoarseSolver="mumps",
    PETSc_options=[
        "-snes_ksp_ew",                     # EW scheme for the KSP
        "-snes_ksp_ew_rtolmax 1e-4",
        "-snes_npicard 2",                  # Number of Picard iterations
        "-snes_monitor",                    # Monitor the nonlinear solver
        "-snes_rtol 1e-5",                  # Relative tolerance for convergence
        "-snes_atol 1e-4",                  # Absolute tolerance for convergence
        "-snes_stol 1e-16",                 # Tolerance for the step length
        "-snes_max_it 100",                  # Maximum number of nonlinear iterations
        "-snes_max_funcs 500000",           # Maximum number of function evaluations
        "-snes_max_linear_solve_fail 10000",# Maximum number of linear solve failures
        "-snes_PicardSwitchToNewton_rtol 5e-2",   # Relative tolerance for switching from Picard to Newton
        "-snes_NewtonSwitchToPicard_it 30",       # Number of Newton iterations before switching to Picard
        "-snes_NewtonSwitchToPicard_rtol 1.1",    # Relative tolerance for switching from Newton to Picard
        "-snes_linesearch_monitor",         # Monitor the line search https://petsc.org/release/manualpages/SNES/SNESNEWTONLS/
        "-snes_linesearch_type cp",
        "-snes_linesearch_maxstep 1",
   

        "-js_ksp_type fgmres",              # Type of Krylov solver for Jacobian solve
        "-js_ksp_max_it 100",                # Maximum number of iterations for Jacobian solve
        "-js_ksp_converged_reason",         # Print the reason for convergence or divergence of the Jacobian solve
        "-js_ksp_monitor",                  # Monitor the Jacobian solve
        "-js_ksp_atol 1e-5",                # Absolute tolerance for convergence of the Jacobian solve
        "-js_ksp_rtol 1e-5",                # Relative tolerance for convergence of the Jacobian solve
    
        "-pstokes mg"
        "-pcmat_type mono",                 # Type of preconditioner for the matrix
        "-pcmat_no_dev_proj",


        "-da_refine_y 1"
    ]
)

	# # Stokes Preconditioner
	# -pstokes mg					# multigrid

	# # Matrix type
	# -pcmat_type mono			# monolithic matrix [coupled MG solver]
	# -pcmat_no_dev_proj
	# -jp_type mg					# multigrid 

	# # Multigrid preconditioner settings
	# -gmg_pc_type mg
	# -gmg_pc_mg_levels 5			# 5 MG levels, gives coarse grid of 32x16x8	
	# -gmg_pc_mg_galerkin			# 	
	# -gmg_pc_mg_type multiplicative
	# -gmg_pc_mg_cycle_type v
	# #-gmg_pc_view 
	# -gmg_pc_mg_log				# monitors time spend in multigrid if using -log_summary @ the end
	
	
	# # RICHARDSON/JACOBI MG Smoothener - [Anton's favorite options]
	# -gmg_mg_levels_ksp_type richardson
	# -gmg_mg_levels_ksp_richardson_scale 0.5
	# -gmg_mg_levels_pc_type jacobi
	# -gmg_mg_levels_ksp_max_it 20
	

	# # DIRECT, REDUNDANT COARSE SOLVER
    #     -crs_ksp_type preonly
    #     -crs_pc_type redundant
	# -crs_pc_redundant_number 8								# split domain in 4 pieces so ever direct solver step is done on 512 cores
	# -crs_redundant_pc_factor_mat_solver_package superlu_dist
# model.Solver=solverMultigrid
model.Solver = sloverDirect 

model.Output = Output(
    # out_dir="model_r1_2.5e15_SA",
    out_surf=1, out_surf_pvd=1,
    # out_mark=1, out_mark_pvd=1,
    out_dev_stress=1, out_strain_rate=1,
    out_surf_topography=1, out_surf_velocity=1, out_surf_amplitude=1)

prepare_lamem(model, 64)#npc
# run_lamem(model, 64)#RUN
