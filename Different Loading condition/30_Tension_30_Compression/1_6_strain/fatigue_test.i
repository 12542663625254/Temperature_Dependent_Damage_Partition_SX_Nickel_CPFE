[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
[]

[Mesh]
  displacements = 'disp_x disp_y disp_z'

  [./gen]
    type = GeneratedMeshGenerator
    dim = 3
    nx = 2
    ny = 2
    nz = 2
    xmin = 0
    xmax = 1
    ymin = 0
    ymax = 1
    zmin = 0
    zmax = 1
  [../]

  [./bot_corner]
    type = ExtraNodesetGenerator
    new_boundary = bot_corner
    input = gen
    coord = '0 0 0'
  [../]

  [./add_side_sets]
    type = SideSetsFromNormalsGenerator
    normals = '1 0 0
               0 1 0
               0 0 1
              -1 0 0
               0 -1 0
               0 0 -1'
    fixed_normal = false
    new_boundary = 'xp_face yp_face zp_face xn_face yn_face zn_face'
    input = bot_corner
  [../]
[]



[Variables]
  [./disp_x]
    order = FIRST
    family = LAGRANGE
    scaling = 1e-4
  [../]
  [./disp_y]
    order = FIRST
    family = LAGRANGE
    scaling = 1e-4
  [../]
  [./disp_z]
    order = FIRST
    family = LAGRANGE
    scaling = 1e-4
  [../]
[]

[AuxVariables]
  [./stress_xx]
    order = FIRST
    family = MONOMIAL
  [../]
  [./stress_yy]
    order = FIRST
    family = MONOMIAL
  [../]
  [./stress_zz]
    order = FIRST
    family = MONOMIAL
  [../]
  [./strain_xx]
    order = FIRST
    family = MONOMIAL
  [../]
  [./strain_yy]
    order = FIRST
    family = MONOMIAL
  [../]
  [./strain_zz]
    order = FIRST
    family = MONOMIAL
  [../]
  [./vonmises]
    order = FIRST
    family = MONOMIAL
  [../]
  [./Ep_eff]
    order = FIRST
    family = MONOMIAL
  [../]
  [./rho_for]
    order = FIRST
    family = MONOMIAL
  [../]
  
  [./phi_1]
    order = FIRST
    family = MONOMIAL
  [../]
  [./Phi]
    order = FIRST
    family = MONOMIAL
  [../]
  [./phi_2]
    order = FIRST
    family = MONOMIAL
  [../]
   [./backstress01]
    order = FIRST
    family = MONOMIAL
  [../]
  [./backstress02]
    order = FIRST
    family = MONOMIAL
  [../]
  [./backstress03]
    order = FIRST
    family = MONOMIAL
  [../]
  [./backstress04]
    order = FIRST
    family = MONOMIAL
  [../]
  [./backstress05]
    order = FIRST
    family = MONOMIAL
  [../]
  [./backstress06]
    order = FIRST
    family = MONOMIAL
  [../]
  [./backstress07]
    order = FIRST
    family = MONOMIAL
  [../]
  [./backstress08]
    order = FIRST
    family = MONOMIAL
  [../]
  [./backstress09]
    order = FIRST
    family = MONOMIAL
  [../]
  [./backstress10]
    order = FIRST
    family = MONOMIAL
  [../]
  [./backstress11]
    order = FIRST
    family = MONOMIAL
  [../]
  [./backstress12]
    order = FIRST
    family = MONOMIAL
  [../]
   [./backstress_avg]
    order = FIRST
    family = MONOMIAL
  [../]
  [./gamma_dot_average]
    order = FIRST
    family = MONOMIAL
  [../]
   [./slip_resistance_avg]
    order = FIRST
    family = MONOMIAL
  [../]
  [./F_norm]
    order = FIRST
    family = MONOMIAL
  [../]
#   [./rho_forest_average]
#    order = FIRST
#    family = MONOMIAL
#  [../]
  
  [./ssd]
    order=FIRST
    family=MONOMIAL
  [../]
  
  [./Prager_term]
    order=FIRST
    family=MONOMIAL
  [../]
  
  [./Dynamic_hardening_term]
    order=FIRST
    family=MONOMIAL
  [../]
  
  [./Static_hardening_term]
    order=FIRST
    family=MONOMIAL
  [../]
  
  [./Prager_term_avg]
    order=FIRST
    family=MONOMIAL
  [../]
  
  [./Dynamic_hardening_term_avg]
    order=FIRST
    family=MONOMIAL
  [../]
  
  [./Static_hardening_term_avg]
    order=FIRST
    family=MONOMIAL
  [../]
  
 
  
  [./Energy]
    order=FIRST
    family=MONOMIAL
  [../]
  
   [./entropy_generation_rate]
    order=FIRST
    family=MONOMIAL
  [../]
  
  
  
  
[]

[Physics/SolidMechanics/QuasiStatic]
  [./all]
    strain = FINITE
    incremental = true
    use_finite_deform_jacobian = true
    volumetric_locking_correction = false
  [../]
[]

[AuxKernels]
  [./stress_xx]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_xx
    index_i = 0
    index_j = 0
  [../]
  [./stress_yy]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_yy
    index_i = 1
    index_j = 1
  [../]
  [./stress_zz]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_zz
    index_i = 2
    index_j = 2
  [../]
  [./strain_zz]
    type = RankTwoAux
    rank_two_tensor = total_strain
    variable = strain_zz
    execute_on = timestep_end
    index_i = 2
    index_j = 2
  [../]
  [./strain_yy]
    type = RankTwoAux
    rank_two_tensor = total_strain
    variable = strain_yy
    execute_on = timestep_end
    index_i = 1
    index_j = 1
  [../]
  [./strain_xx]
    type = RankTwoAux
    rank_two_tensor = total_strain
    variable = strain_xx
    execute_on = timestep_end
    index_i = 0
    index_j = 0
  [../]
  [./vonmises]
    type = RankTwoScalarAux
    rank_two_tensor = stress
    variable = vonmises
    execute_on = timestep_end
    scalar_type = VonMisesStress
  [../]
  [./phi_1]
    type = OutputEuler
    variable = phi_1
    angle_id = 0
    execute_on = timestep_end
  [../]
  [./Phi]
    type = OutputEuler
    variable = Phi
    angle_id = 1
    execute_on = timestep_end
  [../]
  [./phi_2]
    type = OutputEuler
    variable = phi_2
    angle_id = 2
    execute_on = timestep_end
  [../]
  [./Ep_eff]
    type = StateVariable
    variable = Ep_eff
    sdv_id = 47
    execute_on = timestep_end
  [../]
 # [./rho_for]
 #   type = StateVariable
 #   variable = rho_for
 #   sdv_id = 70
 #   execute_on = timestep_end
 # [../]
   [./backstress01]
    type = StateVariable
    variable = backstress01
    sdv_id = 48
    execute_on = timestep_end
  [../]
  [./backstress02]
    type = StateVariable
    variable = backstress02
    sdv_id = 49
    execute_on = timestep_end
  [../]
  [./backstress03]
    type = StateVariable
    variable = backstress03
    sdv_id = 50
    execute_on = timestep_end
  [../]
  [./backstress04]
    type = StateVariable
    variable = backstress04
    sdv_id = 51
    execute_on = timestep_end
  [../]
  [./backstress05]
    type = StateVariable
    variable = backstress05
    sdv_id = 52
    execute_on = timestep_end
  [../]
  [./backstress06]
    type = StateVariable
    variable = backstress06
    sdv_id = 53
    execute_on = timestep_end
  [../]
  [./backstress07]
    type = StateVariable
    variable = backstress07
    sdv_id = 54
    execute_on = timestep_end
  [../]
  [./backstress08]
    type = StateVariable
    variable = backstress08
    sdv_id = 55
    execute_on = timestep_end
  [../]
  [./backstress09]
    type = StateVariable
    variable = backstress09
    sdv_id = 56
    execute_on = timestep_end
  [../]
  [./backstress10]
    type = StateVariable
    variable = backstress10
    sdv_id = 57
    execute_on = timestep_end
  [../]
  [./backstress11]
    type = StateVariable
    variable = backstress11
    sdv_id = 58
    execute_on = timestep_end
  [../]
  [./backstress12]
    type = StateVariable
    variable = backstress12
    sdv_id = 59
    execute_on = timestep_end
  [../]
  
  [./backstress_avg]
    type = StateVariable
    variable = backstress_avg
    sdv_id = 60
    execute_on = timestep_end
  [../]
   [./gamma_dot_average]
    type = StateVariable
    variable = gamma_dot_average
    sdv_id = 61
    execute_on = timestep_end
  [../]
  [./slip_resistance_avg]
    type = StateVariable
    variable = slip_resistance_avg
    sdv_id = 62
    execute_on = timestep_end
  [../]
  [./F_norm]
    type = StateVariable
    variable = F_norm
    sdv_id = 63
    execute_on = timestep_end
  [../]
 #  [./rho_forest_average]
 #   type = StateVariable
 #   variable = rho_forest_average
 #   sdv_id = 77
 #   execute_on = timestep_end
 # [../]
  
  [./ssd]
    type = StateVariable
    variable = ssd
    sdv_id = 76
    execute_on = timestep_end
  [../]
  
  
   [./Prager_term]
    type = StateVariable
    variable = Prager_term
    sdv_id = 78
    execute_on = timestep_end
  [../]
  
   [./Dynamic_hardening_term]
    type = StateVariable
    variable = Dynamic_hardening_term
    sdv_id = 99
    execute_on = timestep_end
  [../]
  
   [./Static_hardening_term]
    type = StateVariable
    variable = Static_hardening_term
    sdv_id = 108
    execute_on = timestep_end
  [../]
  
   [./Prager_term_avg]
    type = StateVariable
    variable = Prager_term_avg
    sdv_id = 113
    execute_on = timestep_end
  [../]
  
   [./Dynamic_hardening_term_avg]
    type = StateVariable
    variable = Dynamic_hardening_term_avg
    sdv_id = 114
    execute_on = timestep_end
  [../]
  
   [./Static_hardening_term_avg]
    type = StateVariable
    variable = Static_hardening_term_avg
    sdv_id = 115
    execute_on = timestep_end
  [../]
  
  [./Energy]
    type = StateVariable
    variable = Energy
    sdv_id = 117
    execute_on = timestep_end
  [../]
  
  
  [./entropy_generation_rate ]
    type = StateVariable
    variable = entropy_generation_rate
    sdv_id = 128
    execute_on = timestep_end
  [../]
  
  
  
  
  
  
[]

[UserObjects]
  [./euler_angle]
    type = EulerAngleReader
    file_name = orientations.in # orientations.in
    execute_on = 'initial'
  [../]
[]

[Functions]

########################## Tension Test for fitting  we use #############################################
# # This function defines a constant velocity for a strain rate of 1e-3/s
# [./top_push]
#   type = PiecewiseLinear
#   x = '0 15  '
#   y = '0 0.015 '
#   [../]
#   [./dts]
#    type = PiecewiseLinear
#    x = '0 0.001'
#    y = '0.00005 0.01'
#  [../]
  
####################### cyclic test (Fatigue) for 1 cycle ######################################
######################## upto 1.5% #################################################
# [./top_push]
# type = PiecewiseLinear
#  x = '0 15  45  75'
#  y = '0 0.015 -0.015  0.015 '
#  [../]
#  [./dts]
#   type = PiecewiseLinear
#   x = '0 0.001'
#   y = '0.00005 0.01'
# [../]
 

####################### cyclic test (Fatigue) for 10 cycles for 0.8% ######################################
[./top_push]
    type = PiecewiseLinear
    
    # Time points (0.8% Amp, 8s Initial Ramp, 30s Holds, 16s Ramps)
    # Cycle Period = 92s (30+16+30+16)
    x = '0 8 38 54 84 100 130 146 176 192 222 238 268 284 314 330 360 376 406 422 452 468 498 514 544 560 590 606 636 652 682 698 728 744 774 790 820 836 866 882 912 928'
    
    # Strain points (0.008 Amplitude)
    y = '0 0.008 0.008 -0.008 -0.008 0.008 0.008 -0.008 -0.008 0.008 0.008 -0.008 -0.008 0.008 0.008 -0.008 -0.008 0.008 0.008 -0.008 -0.008 0.008 0.008 -0.008 -0.008 0.008 0.008 -0.008 -0.008 0.008 0.008 -0.008 -0.008 0.008 0.008 -0.008 -0.008 0.008 0.008 -0.008 -0.008 0.008'
[../]
   [./dts]
  type = PiecewiseLinear
   x = '0 0.001'
   y = '0.00005 0.01'
   [../] 
   
  
######################## Cyclic Test upto 1.0% #################################################
#[./top_push]
#  type = PiecewiseLinear
#  x = '0 10 30 50 '
#  y = '0 0.01 -0.01  0.01 '
#  [../]
#  [./dts]
#  type = PiecewiseLinear
#   x = '0 0.001'
#   y = '0.00005 0.01'
# [../]

####################### 'Dwell in Tension' Cyclic Test 1.2% strain for 60 sec hold time in Tesnion only  #########################################
#[./top_push]
#    type = PiecewiseLinear
#    x = '0 12 72 96 120 180 204 228 288 312 336 396 420 444 504 528 552 612 636 660 720 744 768 828 852 876 936 960 984 1044 1068 1092'
#    y = '0 0.012 0.012 -0.012 0.012 0.012 -0.012 0.012 0.012 -0.012 0.012 0.012 -0.012 0.012 0.012 -0.012 0.012 0.012 -0.012 0.012 0.012 -0.012 0.012 0.012 -0.012 0.012 0.012 -0.012 0.012 0.012 #-0.012 0.012'
#  [../]
#  [./dts]
#  type = PiecewiseLinear
#   x = '0 0.001'
#   y = '0.00005 0.01'
# [../]


####################### 'Dwell in Compression' Cyclic Test 1.2% strain for 60 sec hold time in Tesnion only  #########################################
#[./top_push]
#    type = PiecewiseLinear
#    x = '0 12 36 96 120 144 204 228 252 312 336 360 420 444 468 528 552 576 636 660 684 744 768 792 852 876 900 960 984 1008 1068 1092'
#    y = '0 0.012 -0.012 -0.012 0.012 -0.012 -0.012 0.012 -0.012 -0.012 0.012 -0.012 -0.012 0.012 -0.012 -0.012 0.012 -0.012 -0.012 0.012 -0.012 -0.012 0.012 -0.012 -0.012 0.012 -0.012 -0.012 0.012 #-0.012 -0.012 0.012'
#  [../]
#  [./dts]
#  type = PiecewiseLinear
#   x = '0 0.001'
#   y = '0.00005 0.01'
# [../]


####################### Creep_Fatigue Test 1.2% strain for 30 sec hold time in Tension #########################################

#[./top_push]
#   type = PiecewiseLinear
#    x = '0 12 42 66 90 120 144 168 198 222 246 276 300 324 354 378 402 432 456 480 510 534 558 #588 612 636 666 690 714 744 768 792'
#    y = '0 0.012 0.012 -0.012 0.012 0.012 -0.012 0.012 0.012 -0.012 0.012 0.012 -0.012 0.012 #0.012 -0.012 0.012 0.012 -0.012 0.012 0.012 -0.012 0.012 0.012 -0.012 0.012 0.012 -0.012 0.012 #0.012 
#    -0.012 0.012'
#  [../]
#  [./dts]
#  type = PiecewiseLinear
#   x = '0 0.001'
#   y = '0.00005 0.01'
# [../]

####################### Creep_Fatigue Test 1.2% strain for 30 sec hold time #########################################
#[./top_push]
#    type = PiecewiseLinear
#    x = '0 12 42 66 96 120 150 174 204 228 258 282 312 336 366 390 420 444 474 498 528 552 582 606 636 660 690 714 744 768 798 822 852 876 906 930 960 984 1014 1038 1068 1092'
#    y = '0 0.012 0.012 -0.012 -0.012 0.012 0.012 -0.012 -0.012 0.012 0.012 -0.012 -0.012 0.012 0.012 -0.012 -0.012 0.012 0.012 -0.012 -0.012 0.012 0.012 -0.012 -0.012 0.012 0.012 -0.012 -0.012 #0.012 0.012 -0.012 -0.012 0.012 0.012 -0.012 -0.012 0.012 0.012 -0.012 -0.012 0.012'
#  [../]
#  [./dts]
#  type = PiecewiseLinear
#   x = '0 0.001'
#   y = '0.00005 0.01'
# [../]

####################### Creep_Fatigue Test 1.2% strain for 20 sec hold time #########################################
#[./top_push]
#  type = PiecewiseLinear
#  x = '0  4  24  32  52  60  80   88  108  116  '
#  y = '0 0.012 0.012  -0.012 -0.012  0.012 0.012 -0.012 -0.012 0.012'
#  [../]
#  [./dts]
#  type = PiecewiseLinear
#   x = '0 0.001'
#    y = '0.00005 0.01'
# [../]



#[./top_push]
#  type = PiecewiseLinear
#  x = '0 4 7 14.2 17.2 24.4 27.4 34.6 37.6 44.8 47.8 55 58 65.2 68.2 75.4 78.4 85.6 88.6 95.8 98.8 106 109 116.2 119.2 126.4 129.4 136.6 139.6 146.8 149.8 157 160 167.2 170.2 177.4 180.4 187.6 #190.6 197.8 200.8 208'
#  y = '0 0.012 0.012 -0.012 -0.012 0.012 0.012 -0.012 -0.012 0.012 0.012 -0.012 -0.012 0.012 0.012 -0.012 -0.012 0.012 0.012 -0.012 -0.012 0.012 0.012 -0.012 -0.012 0.012 0.012 -0.012 -0.012 0.012 #0.012 -0.012 -0.012 0.012 0.012 -0.012 -0.012 0.012 0.012 -0.012 -0.012 0.012'
#  [../]
#  [./dts]
#  type = PiecewiseLinear
#   x = '0 0.001'
#    y = '0.00005 0.01'
# [../]


####################### cyclic test (Fatigue) for stress amplitude plot ######################################
######################## upto 1.5% #################################################
#[./top_push]
#  type = PiecewiseLinear
#  x = '0 4 12 20 28 36 44 52 60 68 76 84 92 100 108 116 124 132 140 148 156 164'
#  y = '0 0.012 -0.012 0.012 -0.012 0.012 -0.012 0.012 -0.012 0.012 -0.012 0.012 -0.012 0.012 -0.012 0.012 -0.012 0.012 -0.012 0.012 -0.012 0.012'
# [../]
#  [./dts]
#   type = PiecewiseLinear
#   x = '0 0.001'
#   y = '0.00005 0.01'
# [../]  
#[]
####################### cyclic test (Fatigue) for stress amplitude plot ######################################
######################## upto 1.0% #################################################

#[./top_push]
#    type = PiecewiseLinear
#    x = '0 3.333333 10 16.66666 23.33333 30 36.66666 43.33333 50 56.66666 63.33333 70 76.66666 83.33333 90 96.66666 103.33333 110 116.66666 123.33333 130 136.66666'
#   y = '0 0.011 -0.011 0.011 -0.011 0.011 -0.011 0.011 -0.011 0.011 -0.011 0.011 -0.011 0.011 -0.011 0.011 -0.011 0.011 -0.011 0.011 -0.011 0.011'
#  [../]
#  [./dts]
#    type = PiecewiseLinear
#    x = '0 0.001'
#    y = '0.00005 0.01'
#  [../] 

[]




[BCs]

##############################################
# X : Left and Right
# Y : Bottom and Top
# Z : Back and Front
##############################################
# Plane strain compression
# Compression along z-direction on 'front' boundary
# Constraint along x-direction on 'left' and 'right' boundaries
# y-direction bottom fixed
##############################################

  # --- Fix y at bottom ---
  [./y_roller]
    type = DirichletBC
    variable = disp_y
    boundary = bottom
    value = 0.0
  [../]

  # --- Fix x on left and right (roller BCs) ---
  [./x_roller1]
    type = DirichletBC
    variable = disp_x
    boundary = right
    value = 0.0
  [../]
  [./x_roller2]
    type = DirichletBC
    variable = disp_x
    boundary = left
    value = 0.0
  [../]

  # --- Apply top_push along z (front face) ---
  [./z_push_function]
    type = FunctionDirichletBC
    variable = disp_z
    boundary = front
    function = top_push
  [../]

  # --- Optionally fix back face in z (roller) ---
  [./z_roller]
    type = DirichletBC
    variable = disp_z
    boundary = back
    value = 0.0
  [../]

[]

[Materials]
   [./Pin_Lu_Full_Model]
    type = Pin_Lu_Full_Model
    propsFile = fcc_props.in
    slipSysFile = fcc_slip_sys.in
    num_slip_sys = 12
    num_state_vars = 128 # number of internal state variable
    num_props = 30
    temp = 1033 # K
    tol = 5e-4
    EulerAngFileReader = euler_angle
  [../]
  [./elasticity_tensor]
    type = ComputeCPElasticityTensor
  [../]
[]

[Preconditioning]
  [./SMP]
    type = SMP
    full=true
  [../]
[]

[Executioner]
  type = Transient

  solve_type = 'NEWTON'

  petsc_options = '-snes_ksp_ew'
  petsc_options_iname = '-pc_type -pc_factor_mat_solver_package'
  petsc_options_value = 'lu superlu_dist'

  l_tol = 1e-8
  nl_abs_tol = 5e-7
  nl_rel_tol = 1e-6
  nl_max_its = 20
  nl_forced_its = 1
  l_max_its = 10
  start_time = 0.0
  end_time = 1000000

  [./TimeStepper]
    type = FunctionDT
    function = dts
    min_dt = 1e-8
    cutback_factor_at_failure = 0.1
    growth_factor = 2
  [../]

  [./Predictor]
    type = SimplePredictor
    scale = 1
  [../]

[]

[Postprocessors]
  [./stress_xx]
    type = ElementAverageValue
    variable = stress_xx
  [../]
  [./stress_yy]
    type = ElementAverageValue
    variable = stress_yy
  [../]
  [./stress_zz]
    type = ElementAverageValue
    variable = stress_zz
  [../]
  [./strain_zz]
    type = ElementAverageValue
    variable = strain_zz
  [../]
  [./strain_yy]
    type = ElementAverageValue
    variable = strain_yy
  [../]
  [./strain_xx]
    type = ElementAverageValue
    variable = strain_xx
  [../]
  [./vonmises]
    type = ElementAverageValue
    variable = vonmises
  [../]

  [./phi_1]
    type = ElementAverageValue
    variable = phi_1
  [../]
  [./Phi]
    type = ElementAverageValue
    variable = Phi
  [../]
  [./phi_2]
    type = ElementAverageValue
    variable = phi_2
  [../]
  [./Ep_eff]
    type = ElementAverageValue
    variable = Ep_eff
  [../]
#  [./rho_for]
#    type = ElementAverageValue
#    variable = rho_for
#  [../]
 [./backstress01]
    type = ElementAverageValue
    variable = backstress01
  [../]
  [./backstress02]
    type = ElementAverageValue
    variable = backstress02
  [../]
  [./backstress03]
    type = ElementAverageValue
    variable = backstress03
  [../]
  [./backstress04]
    type = ElementAverageValue
    variable = backstress04
  [../]
  [./backstress05]
    type = ElementAverageValue
    variable = backstress05
  [../]
  [./backstress06]
    type = ElementAverageValue
    variable = backstress06
  [../]
  [./backstress07]
    type = ElementAverageValue
    variable = backstress07
  [../]
  [./backstress08]
    type = ElementAverageValue
    variable = backstress08
  [../]
  [./backstress09]
    type = ElementAverageValue
    variable = backstress09
  [../]
  [./backstress10]
    type = ElementAverageValue
    variable = backstress10
  [../]
  [./backstress11]
    type = ElementAverageValue
    variable = backstress11
  [../]
  [./backstress12]
    type = ElementAverageValue
    variable = backstress12
  [../]
  [./backstress_avg]
    type = ElementAverageValue
    variable = backstress_avg
  [../]
  [./gamma_dot_average]
    type = ElementAverageValue
    variable = gamma_dot_average
  [../]
   [./slip_resistance_avg]
   type = ElementAverageValue
   variable = slip_resistance_avg
 [../]
 [./F_norm]
    type = ElementAverageValue
    variable = F_norm
  [../]
 # [./rho_forest_average]
 #   type = ElementAverageValue
 #   variable = rho_forest_average
 # [../]
   [./ssd]
    type = ElementAverageValue
    variable = ssd
  [../]
  
  
    [./Prager_term]
    type = ElementAverageValue
    variable = Prager_term
  [../]
   
    [./Dynamic_hardening_term]
    type = ElementAverageValue
    variable = Dynamic_hardening_term
  [../]
  
  
    [./Static_hardening_term]
    type = ElementAverageValue
    variable = Static_hardening_term
  [../]
  
    [./Prager_term_avg]
    type = ElementAverageValue
    variable = Prager_term_avg
  [../]
  
   [./Dynamic_hardening_term_avg]
    type = ElementAverageValue
    variable = Dynamic_hardening_term_avg
  [../]
  
   [./Static_hardening_term_avg]
    type = ElementAverageValue
    variable = Static_hardening_term_avg
  [../]
  
   [./Energy]
    type = ElementAverageValue
    variable = Energy
  [../]
   
    [./ entropy_generation_rate]
    type = ElementAverageValue
    variable = entropy_generation_rate
  [../]
 
[]

[Outputs]
  file_base = Fatigue
  csv = true
  print_linear_residuals = true
  perf_graph = true
  time_step_interval = 10
  [./exodus]
   type = Exodus
   time_step_interval = 1
  [../]
  [./cp]
    type = Checkpoint
    time_step_interval = 100
    num_files = 2
  [../]
[]
