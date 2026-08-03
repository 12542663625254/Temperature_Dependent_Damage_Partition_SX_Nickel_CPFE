#include "Pin_Lu_Full_Model.h"
#include "MatrixTools.h"
#include "MooseRandom.h"
#include "MooseException.h"
#include <algorithm>
#include <dlfcn.h>
#include <fstream>
#define QUOTE(macro) stringifyName(macro)

registerMooseObject("RhocpApp", Pin_Lu_Full_Model);registerMooseObject("RhocpApp", Pin_Lu_Full_Model);


InputParameters
Pin_Lu_Full_Model::validParams()
{
  InputParameters params = ComputeStressBase::validParams();
  params.addClassDescription("Stress calculation using the Crystal Plasticity material model for FCC and BCC crystals");
  params.addRequiredParam<FileName>("propsFile", "The file with the material parameters");
  params.addRequiredParam<FileName>("slipSysFile", "The file with the crystallography of slip systems");
  params.addRequiredParam<unsigned int>("num_props", "The number of material properties this UMAT is going to use");
  params.addRequiredParam<unsigned int>("num_slip_sys", "The number of slip systems");
  params.addRequiredParam<unsigned int>("num_state_vars", "The number of state variables this UMAT is going to use");
  params.addParam<Real>("tol", 1.0e-6, "Tolerance");
  params.addCoupledVar("temp", 300, "Temperature");
  params.addParam<UserObjectName>("EulerAngFileReader", "Name of the EulerAngleReader UO");
  params.addParam<UserObjectName>("EBSDFileReader", "Name of the EBSDReader UO");
  params.addParam<UserObjectName>("GrainAreaSize", "Name of the GrainAreaSize UO");
  params.addParam<int>("isEulerRadian", 0, "Are Euler angles specified in radians");
  params.addParam<int>("isEulerBunge", 0, "Are Euler angles specified in Bunge notation");
  params.addParam<std::vector<MaterialPropertyName>>(
      "eigenstrain_names", {}, "List of eigenstrains to be applied in this strain calculation");
  return params;
}

Pin_Lu_Full_Model::Pin_Lu_Full_Model(const InputParameters & parameters) :
    ComputeStressBase(parameters),
    _propsFile(getParam<FileName>("propsFile")),
    _slipSysFile(getParam<FileName>("slipSysFile")),
    _num_props(getParam<unsigned int>("num_props")),
    _num_slip_sys(getParam<unsigned int>("num_slip_sys")),
    _num_state_vars(getParam<unsigned int>("num_state_vars")),
    _tol(getParam<Real>("tol")),
    _temp(coupledValue("temp")),
    _EulerAngFileReader(isParamValid("EulerAngFileReader")
                               ? &getUserObject<EulerAngleReader>("EulerAngFileReader")
                               : NULL),
    _EBSDFileReader(isParamValid("EBSDFileReader")
                               ? &getUserObject<EBSDMeshReader>("EBSDFileReader")
                               : NULL),
    _GrainAreaSize(isParamValid("GrainAreaSize")
                               ? &getUserObject<GrainAreaSize>("GrainAreaSize")
                               : NULL),
    _isEulerRadian(getParam<int>("isEulerRadian")),
    _isEulerBunge(getParam<int>("isEulerBunge")),
    _euler_ang(declareProperty<Point>("euler_ang")),
    _eigenstrain_names(getParam<std::vector<MaterialPropertyName>>("eigenstrain_names")),
    _eigenstrains(_eigenstrain_names.size()),
    _eigenstrains_old(_eigenstrain_names.size()),
    _state_var(declareProperty<std::vector<Real> >("state_var")),
    _state_var_old(getMaterialPropertyOld<std::vector<Real> >("state_var")),
    _properties(declareProperty<std::vector<Real> >("properties")),
    _properties_old(getMaterialPropertyOld<std::vector<Real> >("properties")),
    _y(declareProperty<std::vector<std::vector<Real> > >("slip plane normals")),
    _y_old(getMaterialPropertyOld<std::vector<std::vector<Real> > >("slip plane normals")),
    _z(declareProperty<std::vector<std::vector<Real> > >("slip plane directions")),
    _z_old(getMaterialPropertyOld<std::vector<std::vector<Real> > >("slip plane directions")),
    _deformation_gradient(getMaterialProperty<RankTwoTensor>("deformation_gradient")),
    _deformation_gradient_old(getMaterialPropertyOld<RankTwoTensor>("deformation_gradient")),
    _strain_increment(getMaterialPropertyByName<RankTwoTensor>(_base_name + "strain_increment")),
    _rotation_increment(getMaterialPropertyByName<RankTwoTensor>(_base_name + "rotation_increment")),
    _stress_old(getMaterialPropertyOld<RankTwoTensor>(_base_name + "stress")),
    _Cel_cp(declareProperty<RankFourTensor>("Cel_cp"))
{
  for (auto i : make_range(_eigenstrain_names.size())) {
    _eigenstrains[i] = &getMaterialProperty<RankTwoTensor>(_eigenstrain_names[i]);
    _eigenstrains_old[i] = &getMaterialPropertyOld<RankTwoTensor>(_eigenstrain_names[i]);
  }
}

Pin_Lu_Full_Model::~Pin_Lu_Full_Model()
{
}

void Pin_Lu_Full_Model::initQpStatefulProperties()
{
  ComputeStressBase::initQpStatefulProperties();

  //Initialize state variable vector
  _state_var[_qp].resize(_num_state_vars);

  // Initialize _properties vector
  _properties[_qp].resize(_num_props);

  // Initialize _y and _z vectors
  _y[_qp].resize(3, std::vector<Real>(_num_slip_sys));
  _z[_qp].resize(3, std::vector<Real>(_num_slip_sys));

} // End initQpStatefulProperties()

void Pin_Lu_Full_Model::computeQpStress()
{
  Real PI = 4.0*atan(1.0);
  Real psi[3]; // Euler angles
  RankTwoTensor E_el; // Elastic Green strain tensor
  RankTwoTensor E_tot; // Green strain tensor E_el + E_p
  RankTwoTensor F0; // F at beginning of sub increment
  RankTwoTensor F1; // F at end of sub increment
  RankTwoTensor F_el; //  Elastic part of F
  RankTwoTensor F_el_inv; // Inverse of elastic part of F
  RankTwoTensor F_p_inv_0; // Inverse of plastic part of F at beginning of step
  RankTwoTensor F_p_inv; // Inverse of plastic part of F at end of step
  RankTwoTensor E_p; // Plastic strain tensor
  Real backstress_avg;       // average backstress value
  Real rho_forest_average;   // average forest density
  Real gamma_dot_average;        // average gamma_dot value
  Real slip_resistance_avg;  // average slip resistance value
  Real slip_resistance_avg0; //
  Real F_norm;               // Modulus F at end of sub increment
  Real E_eff;
  Real E_p_eff;
  Real E_p_eff_cum;
  

  // entropy generation rate term
 Real entropy_generation_rate_0,entropy_generation_rate;
 



 // Real rho_for_avg;
  Real rho_SSD_avg0, rho_SSD_avg;

 //Backstress term's average value
  Real Prager_term_avg_0,Prager_term_avg;
  Real Dynamic_hardening_term_avg_0,Dynamic_hardening_term_avg;
  Real Static_hardening_term_avg_0,Static_hardening_term_avg;
  // KM equation terms
  Real First_term_KM_avg_0,First_term_KM_avg;
  Real Second_term_KM_avg_0,Second_term_KM_avg;




  Real dir_cos[3][3]; // Rotation matrix
  Real dir_cos0[3][3]; // Rotation matrix at the beginning of simulation
  Real C0[3][3][3][3]; // Elastic stiffness tensor
  Real C[3][3][3][3];
  Real C_avg[3][3][3][3];

  RankTwoTensor sig_avg;  //  Average Stress 
  RankTwoTensor sig;      // Cauchy stress
  RankTwoTensor Spk2;     // second piola kirchoff stress

  Real ddpdsig[3][3][3][3];
  Real ddsdde_4th[3][3][3][3];

  // slip system dependent arrays
  std::vector<std::vector<Real>> xs0(3, std::vector<Real>(_num_slip_sys)); // intermediate config slip directions in global coords
  std::vector<std::vector<Real>> xs(3, std::vector<Real>(_num_slip_sys)); // current config slip directions in global coords
  std::vector<std::vector<Real>> xm0(3, std::vector<Real>(_num_slip_sys)); // intermediate config plane normals in global coords
  std::vector<std::vector<Real>> xm(3, std::vector<Real>(_num_slip_sys)); // current config plane normals in global coords

  std::vector<Real> rho_SSD0(_num_slip_sys); // SSD density at beginning of step
  std::vector<Real> rho_SSD(_num_slip_sys); // SSD density at end of step

  std::vector<Real> rho_for_0(_num_slip_sys); // mobile dislocation density at beginning of step ...... not usedin this model
  std::vector<Real> rho_for(_num_slip_sys); // mobile dislocation density at end of step ........ not usedin this model
  std::vector<Real> bstress0(_num_slip_sys); // backstress at beginning of step
  std::vector<Real> bstress(_num_slip_sys); // backstress at end of step
  std::vector<Real> tau(_num_slip_sys); // resolved shear stress
  std::vector<Real> tau_eff(_num_slip_sys); // effective resolved shear stress
  std::vector<Real> s_a(_num_slip_sys); // athermal slip resistance
  std::vector<Real> s_a_0(_num_slip_sys); // athermal slip resistance
  std::vector<Real> s_t(_num_slip_sys); // thermal slip resistance
  std::vector<Real> gamma_dot(_num_slip_sys); // shear strain rate
  std::vector<Real> gamma_try(_num_slip_sys); // trial shear strain rate
  //std::vector<Real> gamma_try_g(_num_slip_sys); // trial shear strain rate due to glide

  ////////////store the Backstressterm separately///////////////////////////
  
  std::vector<Real> Prager_term_0(_num_slip_sys);  //Intial value of 'First term of the Kinematic hardening equation (Modified A-F backstress equation)'
  std::vector<Real> Prager_term(_num_slip_sys);  // Value of the first term of the Modified A-F backstress
  std::vector<Real> Dynamic_hardening_term_0(_num_slip_sys);  // Initial of the second term of the Modified A-F backstress
  std::vector<Real> Dynamic_hardening_term(_num_slip_sys);  // Value of the second term of the Modified A-F backstress
  std::vector<Real> Static_hardening_term_0(_num_slip_sys);  // Initial of the third term of the Modified A-F backstress
  std::vector<Real> Static_hardening_term(_num_slip_sys);  // Value of the third term of the Modified A-F backstress
  
  // Kocks Mecking terms print
  std::vector<Real> First_term_KM_0(_num_slip_sys);
  std::vector<Real> First_term_KM(_num_slip_sys);
  
  std::vector<Real> Second_term_KM_0(_num_slip_sys);
  std::vector<Real> Second_term_KM(_num_slip_sys);








  // energy by  
  std::vector<Real> Energy_0(_num_slip_sys);  //Intial value of energy 
  std::vector<Real> Energy(_num_slip_sys);  // value of energy



  std::vector<Real> residual(_num_slip_sys); // residual during Newton-Raphson iteration

  std::vector<std::vector<Real>> dTaudGd(_num_slip_sys, std::vector<Real>(_num_slip_sys)); // d(tau)/d(gamma_dot)

  // interaction coefficient martix between different slip systems
  std::vector<std::vector<Real>> A(_num_slip_sys, std::vector<Real>(_num_slip_sys));
  std::vector<std::vector<Real>> H(_num_slip_sys, std::vector<Real>(_num_slip_sys));

  // Our model is based on V.R decomposition as opposed to MOOSE's R.U decomposition of the deformation gradient
  // F = V.R = R.U
  // V = R.U.R_T
  // Our model needs everything in the current configuration
  RankTwoTensor delV;
  delV = _rotation_increment[_qp] * _strain_increment[_qp] * _rotation_increment[_qp].transpose();

  // Rotate stress to current configuration
  _stress[_qp] = _rotation_increment[_qp] * _stress_old[_qp] * _rotation_increment[_qp].transpose();

  // Recover "old" state variables
  std::copy(_state_var_old[_qp].begin(), _state_var_old[_qp].end(), _state_var[_qp].begin());
  std::copy(_properties_old[_qp].begin(), _properties_old[_qp].end(), _properties[_qp].begin());

  // Recover "old" miller indices
  std::copy(_y_old[_qp].begin(), _y_old[_qp].end(), _y[_qp].begin());
  std::copy(_z_old[_qp].begin(), _z_old[_qp].end(), _z[_qp].begin());

  if (_t_step <= 1) {
  // read in Euler angles
  if (_EBSDFileReader){
      Point p = _current_elem->vertex_average();
      EBSDAccessFunctors::EBSDPointData data = _EBSDFileReader->getData(p);
      _grainid = data._feature_id;

      psi[0] = data._phi1;
      psi[1] = data._Phi;
      psi[2] = data._phi2;
  }
  else if(_EulerAngFileReader){
      _grainid = _current_elem->subdomain_id();

      psi[0] = _EulerAngFileReader->getData(_grainid,0);
      psi[1] = _EulerAngFileReader->getData(_grainid,1);
      psi[2] = _EulerAngFileReader->getData(_grainid,2);
  }
  else{
      mooseError("Euler angle data not found.");
  }

  // read grain size, available
  if (_GrainAreaSize){
    _grain_size = _GrainAreaSize->getGrainSize(_grainid);
  }
  else{
    _grain_size = 1.e10;
  }
  }

  // Read in material parameters
  if (_t_step <= 1){
      readPropsFile();
  }
  assignProperties();

  // Initialize Kronecker delta tensor
  Real del[3][3];
  for (unsigned int i = 0; i < 3; i++){
    for (unsigned int j = 0; j < 3; j++){
      del[i][j] = 0;
    }
    del[i][i] = 1;
  }

  // Read in slip system data
  if (_t_step <= 1){
    MooseUtils::checkFileReadable(_slipSysFile);
    std::ifstream file_slip_sys;
    file_slip_sys.open(_slipSysFile.c_str());

    // Assign slip system normals and slip directions
    for (unsigned int i = 0; i < _num_slip_sys; i++) {
      for(unsigned int j = 0; j < 3; j++){
        //y[j][i] = reader.getData(c++)[0];
        file_slip_sys >> _y[_qp][j][i];
      }
      for (unsigned int j = 0; j < 3; j++){
        //z[j][i] = reader.getData(c++)[0];
        file_slip_sys >> _z[_qp][j][i];
      }
    }
    file_slip_sys.close();

    // Normalize the Miller indices to unit vectors
    for (unsigned int i = 0; i < _num_slip_sys; i++) {
      normalize_vector(&_y[_qp][0][i],&_y[_qp][1][i],&_y[_qp][2][i]);
      normalize_vector(&_z[_qp][0][i],&_z[_qp][1][i],&_z[_qp][2][i]);
    }

    for (unsigned int k = 0; k < _num_slip_sys; k++) {
      Real sum1 = 0.0;
      for (unsigned int i = 0; i<3; i++) {
        sum1 = sum1 + _y[_qp][i][k]*_z[_qp][i][k];
      }
      if (abs(sum1) > 1.0e-8) {
        mooseError("The Miller indices are WRONG for slip system:",(k + 1));
        break;
      }
    }
  }
   
  // // Initialize interaction coefficient matrix
  // for (unsigned int i = 0; i < _num_slip_sys; i++) {
  //   for (unsigned int j = 0; j < _num_slip_sys; j++) {
  //     A[i][j] = Alatent;
  //     H[i][j] = 0.0;
  //   }
  //   A[i][i] = 1.0;
  //   H[i][i] = 1.0;
  //}
    static const double A_table[12][12] = {
      {0.100000, 0.220000, 0.220000, 0.150627, 0.424262, 0.358637, 0.150627, 0.358637, 0.424262, 0.141421, 0.358637, 0.358637},
      {0.220000, 0.100000, 0.220000, 0.424262, 0.150627, 0.358637, 0.358637, 0.141421, 0.358637, 0.358637, 0.150627, 0.424262},
      {0.220000, 0.220000, 0.100000, 0.358637, 0.358637, 0.141421, 0.424262, 0.358637, 0.150627, 0.358637, 0.424262, 0.150627},
      {0.150627, 0.424262, 0.358637, 0.100000, 0.220000, 0.220000, 0.141421, 0.358637, 0.358637, 0.150627, 0.358637, 0.424262},
      {0.424262, 0.150627, 0.358637, 0.220000, 0.100000, 0.220000, 0.358637, 0.150627, 0.424262, 0.358637, 0.141421, 0.358637},
      {0.358637, 0.358637, 0.141421, 0.220000, 0.220000, 0.100000, 0.358637, 0.424262, 0.150627, 0.424262, 0.358637, 0.150627},
      {0.150627, 0.358637, 0.424262, 0.141421, 0.358637, 0.358637, 0.100000, 0.220000, 0.220000, 0.150627, 0.424262, 0.358637},
      {0.358637, 0.141421, 0.358637, 0.358637, 0.150627, 0.424262, 0.220000, 0.100000, 0.220000, 0.424262, 0.150627, 0.358637},
      {0.424262, 0.358637, 0.150627, 0.358637, 0.424262, 0.150627, 0.220000, 0.220000, 0.100000, 0.358637, 0.358637, 0.141421},
      {0.141421, 0.358637, 0.358637, 0.150627, 0.358637, 0.424262, 0.150627, 0.424262, 0.358637, 0.100000, 0.220000, 0.220000},
      {0.358637, 0.150627, 0.424262, 0.358637, 0.141421, 0.358637, 0.424262, 0.150627, 0.358637, 0.220000, 0.100000, 0.220000},
      {0.358637, 0.424262, 0.150627, 0.424262, 0.358637, 0.150627, 0.358637, 0.358637, 0.141421, 0.220000, 0.220000, 0.100000}};

     // copy element-by-element:
      for (unsigned i = 0; i < _num_slip_sys; ++i)
     {
      for (unsigned j = 0; j < _num_slip_sys; ++j)
     {
      A[i][j] = A_table[i][j]; // scalar-to-scalar
      H[i][j] = 0.0;
     }
      H[i][i] = 1.0;
     }


  if (_t_step == 1) {
    // Convert Euler angles to radians
    psi[0] = psi[0]*PI/180;
    psi[1] = psi[1]*PI/180;
    psi[2] = psi[2]*PI/180;

    // Initialize F_p_inv_0
    for (unsigned int i = 0; i<3; i++) {
      for(unsigned int j = 0; j<3; j++) {
        F_p_inv_0(i,j) = 0.0;
      }
      F_p_inv_0(i,i) = 1.0;
    }

    // Initialize E_p.
    for (unsigned int i = 0; i < 3; i++) {
      for(unsigned int j = 0; j < 3; j++) {
        E_p(i,j) = 0;
      }
    }

    // Initialize E_eff
    E_eff = 0;

    // Initialize E_p_eff
    E_p_eff = 0;

    // Initialize E_p_eff_cum
    E_p_eff_cum = 0;

    // Initialize average values of defect densities and the corresponding loop sizes
   // rho_for_avg0 = rho_for_zero;
    backstress_avg=0.0;       
    gamma_dot_average = 0.0;
    slip_resistance_avg0=0.0;
    F_norm=0.0;
  
    rho_SSD_avg0 = rho_SSD_zero;
    Prager_term_avg_0 = 0;
    Dynamic_hardening_term_avg_0 = 0;
    Static_hardening_term_avg_0 = 0;
    entropy_generation_rate_0=0;
    First_term_KM_avg_0=0;
    Second_term_KM_avg_0=0;

   for (unsigned int n = 0; n<_num_slip_sys; n++){
      rho_SSD0[n] = rho_SSD_zero;
    }
  
    for (unsigned int n = 0; n<_num_slip_sys; n++){
      bstress0[n] = 0.0;
    }
    
    for (unsigned int n = 0; n<_num_slip_sys; n++){
      Prager_term_0[n] = 0.0;
    }

    
    for (unsigned int n = 0; n<_num_slip_sys; n++){
      Dynamic_hardening_term_0[n] = 0.0;
    }

   for (unsigned int n = 0; n<_num_slip_sys; n++){
      Static_hardening_term_0[n] = 0.0;
    }

    for (unsigned int n = 0; n<_num_slip_sys; n++){
      First_term_KM_0[n] = 0.0;
    }

    for (unsigned int n = 0; n<_num_slip_sys; n++){
      Second_term_KM_0[n] = 0.0;
    }




     for (unsigned int n = 0; n<_num_slip_sys; n++){
      Energy_0[n] = 0.0;
    }


  } // time  loop

  //  End of initializations.  Read in internal variables

  if(_t_step > 1) {
    // checkpoint("Initializing ISVs, nonzero time step")
    int n = -1;
    // Read in dir_cos0 SDV 1-9
    for (unsigned int i = 0; i<3; i++) {
      for (unsigned int j = 0; j < 3; j++) {
        n = n + 1;
        dir_cos0[i][j] = _state_var[_qp][n];
      }
    }

    // Read in dir_cos SDV 10-18
    for (unsigned int i = 0; i<3; i++) {
      for (unsigned int j = 0; j < 3; j++) {
        n = n + 1;
        dir_cos[i][j] = _state_var[_qp][n];
      }
    }

    // Read the elastic part of F SDV 19-27
    for (unsigned int i = 0; i < 3; i++) {
      for (unsigned int j = 0; j < 3; j++) {
        n = n + 1;
        // F_el[i][j] = _state_var[_qp][n]; // Not needed. This is only for post-processing
      }
    }

    // Read inverse of the plastic part of F SDV 28-36
    for (unsigned int i = 0; i < 3; i++) {
      for (unsigned int j = 0; j < 3; j++) {
        n = n + 1;
        F_p_inv_0(i,j) = _state_var[_qp][n];
      }
    }

    // Read E_p SDV 37-45
    for (unsigned int i = 0; i < 3; i++) {
      for (unsigned int j = 0; j < 3; j++) {
        n = n + 1;
        E_p(i,j) = _state_var[_qp][n];
      }
    }

    // Read E_eff SDV 46
    n = n + 1;
    E_eff = _state_var[_qp][n];

    // Read E_p_eff SDV 47
    n = n + 1;
    E_p_eff = _state_var[_qp][n];

    
    // Read back stress values SDV 48-59
    for (unsigned int i = 0; i < _num_slip_sys; i++) {
      n = n + 1;
      bstress0[i] = _state_var[_qp][n];
    }
   
    //Read backstress_avg SDV 60
    n=n+1;
    backstress_avg = _state_var[_qp][n];
   
    //Read gamma_dot_avg SDV 61
    n=n+1;
    gamma_dot_average = _state_var[_qp][n];
    // Read Slip_resistance average value; ISV 62
    n=n+1;   
    slip_resistance_avg0=_state_var[_qp][n];
    //Read modulus of deformation gradient SDV 63
    n=n+1;
    F_norm = _state_var[_qp][n];
    
     // Read SSD density values SDV 64-75
    for (unsigned int i = 0; i < _num_slip_sys; i++) {
      n = n + 1;
      rho_SSD0[i] = _state_var[_qp][n];
    }
   
     // Read average values of defect densities and loop sizes SDV 76
    n = n + 1;
    rho_SSD_avg0 = _state_var[_qp][n];
    
    // Read first term of back stress  SDV 77-88
    for (unsigned int i = 0; i < _num_slip_sys; i++) {
      n = n + 1;
      Prager_term_0[i] = _state_var[_qp][n];
    } 
    
     // Read second term of back stress  SDV 89-100
    for (unsigned int i = 0; i < _num_slip_sys; i++) {
      n = n + 1;
      Dynamic_hardening_term_0[i] = _state_var[_qp][n];
    } 
    
    // Read third term of back stress  SDV 101-112
    for (unsigned int i = 0; i < _num_slip_sys; i++) {
      n = n + 1;
      Static_hardening_term_0[i] = _state_var[_qp][n];
    }  

      // Read average values of backstress term SDV 113
    n = n + 1;
    Prager_term_avg_0 = _state_var[_qp][n];

      // Read average values of backstress term SDV 114
    n = n + 1;
    Dynamic_hardening_term_avg_0 = _state_var[_qp][n];
    
     // Read average values of backstress term SDV 115
    n = n + 1;
    Static_hardening_term_avg_0 = _state_var[_qp][n];

     // energy
 
     // Read energy  SDV 116-127
    for (unsigned int i = 0; i < _num_slip_sys; i++) {
      n = n + 1;
      Energy_0[i] = _state_var[_qp][n];
    } 
    
    // entropy generation rate  SVD 128
       
     n=n+1;
     entropy_generation_rate_0=_state_var[_qp][n];

    // Kocks mecking terms ISV  129-140

     for (unsigned int i = 0; i < _num_slip_sys; i++) {
      n = n + 1;
      First_term_KM_0[i] = _state_var[_qp][n];
    } 

    // Kocks mecking terms ISV  141-152
  

     for (unsigned int i = 0; i < _num_slip_sys; i++) {
      n = n + 1;
      Second_term_KM_0[i] = _state_var[_qp][n];
    } 


   // Read average values of KM term SDV 153
    n = n + 1;
    First_term_KM_avg_0 = _state_var[_qp][n];
  
  // Read average values of KM term SDV 154
    n = n + 1;
    Second_term_KM_avg_0 = _state_var[_qp][n];









    // std::cout << "\nTotal no. of state vars read:" << n;
  } // End of initializations

  if(_t_step == 1) {
    Real s1 = sin(psi[0]);
    Real c1 = cos(psi[0]);
    Real s2 = sin(psi[1]);
    Real c2 = cos(psi[1]);
    Real s3 = sin(psi[2]);
    Real c3 = cos(psi[2]);

    dir_cos0[0][0] = c1*c3-s1*s3*c2;
    dir_cos0[0][1] = s1*c3+c1*s3*c2;
    dir_cos0[0][2] = s3*s2;
    dir_cos0[1][0] = -c1*s3-s1*c3*c2;
    dir_cos0[1][1] = -s1*s3+c1*c3*c2;
    dir_cos0[1][2] = c3*s2;
    dir_cos0[2][0] = s1*s2;
    dir_cos0[2][1] = -c1*s2;
    dir_cos0[2][2] = c2;

    for (unsigned int i = 0; i < 3; i++){
      for (unsigned int j = 0; j < 3; j++){
        dir_cos[i][j] = dir_cos0[i][j];
      }
    }
  }

  // Initialize ANISOTROPIC 4th rank elastic stiffness tensor
  for (unsigned int i = 0; i < 3; i++) {
    for (unsigned int j = 0; j < 3; j++) {
      for (unsigned int k = 0; k < 3; k++) {
        for (unsigned int l = 0; l < 3; l++) {
          C0[i][j][k][l] = C12*del[i][j]*del[k][l] + C44*(del[i][k]*del[j][l]+del[i][l]*del[k][j]);
        }
      }
    }
  }
  C0[0][0][0][0] = C11;
  C0[1][1][1][1] = C11;
  C0[2][2][2][2] = C11;

  // Initialize arrays for averaging over grains
  for (unsigned int i = 0; i < 3; i++) {
    for (unsigned int j = 0; j < 3; j++) {
      sig_avg(i,j) = 0.0;
    }
  }

  // Begin calculations over the element
  rotate_4th(dir_cos, C0, C);

  // Convert Miller Indices to global coordinates
  for (unsigned int n = 0; n < _num_slip_sys; n++) {
    for (unsigned int i = 0; i < 3; i++) {
      xs0[i][n] = 0.0;
      xm0[i][n] = 0.0;
      for (unsigned int j = 0; j < 3; j++) {
        xs0[i][n] = xs0[i][n] + dir_cos0[i][j]*_z[_qp][j][n];
        xm0[i][n] = xm0[i][n] + dir_cos0[i][j]*_y[_qp][j][n];
      }
    }
  }

  // Initialize number of sub-increments.  Note that the subincrement initially equals the total increment. This remains the case unless the process starts to diverge.
  unsigned int N_incr, N_incr_total, iNR, N_ctr;
  Real dt_incr;
  RankTwoTensor dfgrd0, dfgrd1, Ftheta, Ftheta_old;
  RankTwoTensor total_eigenstrain, total_eigenstrain_old;
  RankTwoTensor array1, array2;
  Real array1_matrix[3][3], array2_matrix[3][3];

  bool converged, improved;

  //// Initialize deformation gradients for beginning and end of subincrement.
  //for (auto i : make_range(_eigenstrain_names.size())) {
  //  total_eigenstrain += (*_eigenstrains[i])[_qp];
  //  total_eigenstrain_old += (*_eigenstrains_old[i])[_qp];
  //}
  //
  //Ftheta.zero();
  //Ftheta_old.zero();
  //for (int i=0; i<3; ++i)
  //{
  //    Ftheta(i,i) = sqrt(1.0 + 2.0*total_eigenstrain(i,i));
  //    Ftheta_old(i,i) = sqrt(1.0 + 2.0*total_eigenstrain_old(i,i));
  //}
  //
  dfgrd0 = _deformation_gradient_old[_qp];
  dfgrd1 = _deformation_gradient[_qp] ;
  F0 = dfgrd0;
  F1 = dfgrd1;

  // counters for time step sub-increment
  N_incr = 1;
  N_incr_total = 1;

  Real sse_old, sse_ref;

  // Begin time step subincrement
  do {
  dt_incr = _dt/N_incr_total;

  F0 = dfgrd0 + (dfgrd1 - dfgrd0)*(N_incr - 1)/N_incr_total;
  F1 = dfgrd0 + (dfgrd1 - dfgrd0)*N_incr/N_incr_total;

  //Initial vale of slip resistance_avg;
  slip_resistance_avg=slip_resistance_avg0;

 // Initialize mobile dislocation density
  for (unsigned int n = 0; n < _num_slip_sys; n++) {
    rho_SSD[n] = rho_SSD0[n];
  }
 
  // Initialize backstress
  for (unsigned int n = 0; n < _num_slip_sys; n++) {
    bstress[n] = bstress0[n];
  }
 
  // Initialize backstress first term
  for (unsigned int n = 0; n < _num_slip_sys; n++) {
    Prager_term[n] =Prager_term_0[n];
  }
  
   // Initialize backstress second term
  for (unsigned int n = 0; n < _num_slip_sys; n++) {
    Dynamic_hardening_term[n] =Dynamic_hardening_term_0[n];
  }

   // Initialize backstress third term
  for (unsigned int n = 0; n < _num_slip_sys; n++) {
    Static_hardening_term[n] =Static_hardening_term_0[n];
  }

// Initialize energy term
  for (unsigned int n = 0; n < _num_slip_sys; n++) {
    Energy[n] =Energy_0[n];
  }
  

  // Initialize KM terms first term
  for (unsigned int n = 0; n < _num_slip_sys; n++) {
    First_term_KM[n] =First_term_KM_0[n];
  }

   // Initialize KM terms second term
  for (unsigned int n = 0; n < _num_slip_sys; n++) {
    Second_term_KM[n] =Second_term_KM_0[n];
  }




  // Multiply F() by F_p_inv() to get F_el()
  F_el = F0 * F_p_inv_0;
  if (F_el.det() == 0.0e0) {
    F_el(0,0) = 1.0e0;
    F_el(1,1) = 1.0e0;
    F_el(2,2) = 1.0e0;
  }

  F_el_inv = F_el.inverse();

  // Rotate xs0 and xm0 to current coordinates, called xs and xm
  for (unsigned int n = 0; n < _num_slip_sys; n++) {
    for (unsigned int i = 0; i < 3; i++) {
      xs[i][n] = 0.0;
      xm[i][n] = 0.0;
      for (unsigned int j = 0; j < 3; j++) {
        xs[i][n] = xs[i][n] + F_el(i,j)*xs0[j][n];
        xm[i][n] = xm[i][n] + xm0[j][n]*F_el_inv(j,i);
      }
    }
  }

  // Calculate elastic Green strain
  E_el = F_el.transpose() * F_el;
  E_el(0,0) = E_el(0,0) - 1.0;
  E_el(1,1) = E_el(1,1) - 1.0;
  E_el(2,2) = E_el(2,2) - 1.0;
  E_el = 0.5 * E_el;

  // Multiply the anisotropic stiffness tensor by the Green strain to get the 2nd Piola Kirkhhoff stress
  Real E_el_matrix[3][3], Spk2_matrix[3][3];
  for (unsigned int i = 0; i < 3; i++) {
    for (unsigned int j = 0; j < 3; j++) {
      E_el_matrix[i][j] = E_el(i,j);
    }
  }
  aaaa_dot_dot_bb(C,E_el_matrix,Spk2_matrix);
  for (unsigned int i = 0; i < 3; i++) {
    for (unsigned int j = 0; j < 3; j++) {
      Spk2(i,j) = Spk2_matrix[i][j];
    }
  }

  // Convert from PK2 stress to Cauchy stress
  Real det = F_el.det();
  sig =  F_el * Spk2 * F_el.transpose();
  for (unsigned int i = 0; i < 3; i++) {
    for (unsigned int j = 0; j < 3; j++) {
      sig(i,j) = sig(i,j)/det;
    }
  }

  // Calculate resolved shear stress for each slip system
  for (unsigned int k = 0; k < _num_slip_sys; k++) {
    tau[k] = 0.0;
    for (unsigned int j = 0; j < 3; j++) {
      for (unsigned int i = 0; i < 3; i++) {
        tau[k] = tau[k] + xs[i][k]*xm[j][k]*sig(i,j);
      }
    }
  }

  // Calculate effective shear stress for each slip system
  for (unsigned int k = 0; k < _num_slip_sys; k++) {
    tau_eff[k] = tau[k] - bstress[k];
  }

  // Calculate athermal slip resistance for each slip system for Pin lu Model
  for (unsigned int ia = 0; ia < _num_slip_sys; ia++) {
    Real sum1 = 0.0;
    for (unsigned int ib = 0; ib < _num_slip_sys; ib++) {
      sum1 = sum1 + A[ia][ib]*(rho_SSD[ib]);
    }
    s_a[ia] = Slip_resistance_const +k_0*G*b_berg*sqrt(sum1);
  }


   //NOt used here used in Li wang flow rule
  // Calculate reference shear stress for each slip system.
  for (unsigned int ia = 0; ia < _num_slip_sys; ia++) {
    s_t[ia] = frictional_stress;
  }

  // Calculate 1st estimate of gamma_dot for each slip system.

  for(unsigned int k = 0; k < _num_slip_sys; k++) {
    
   // if(abs(tau_eff[k]) > s_a[k]) {
      // flowrule for Li Wang Model 
      //gamma_dot[k] = gammadot0*exp((-delF0/B_k/temperature)*power(1 - power((abs(tau_eff[k]) - s_a[k])/s_t[k], p), q))*sgn(tau_eff[k]);
       
      // Flow rule for Pin Lu Updated Model
      gamma_dot[k] = gammadot0 * power((std::abs(tau_eff[k] / s_a[k])), (1 / m)) * sgn(tau_eff[k]);
     
   //    } else
   //   {
     // gamma_dot[k] = 1e-18;
   //   }
  }

 //************Calculate d(Tau)/d(Gamma_dot), Eqn.3.47 Mcginty******************


  // Calculate d(Tau)/d(Gamma_dot)
  for(unsigned int ia = 0; ia < _num_slip_sys; ia++) {
    for(unsigned int ib = 0; ib < _num_slip_sys; ib++) {
      dTaudGd[ia][ib] = 0.0;

      for(unsigned int i = 0; i < 3; i++) {
        for(unsigned int j = 0; j < 3; j++) {
          array1_matrix[i][j] = xs0[i][ib]*xm0[j][ib];
        }
      }
      aaaa_dot_dot_bb(C0,array1_matrix,array2_matrix);

      for(unsigned int i = 0; i < 3; i++) {
        for(unsigned int j = 0; j < 3; j++) {
          array1_matrix[i][j] = xs0[i][ia]*xm0[j][ia];
        }
      }
      dTaudGd[ia][ib] = aa_dot_dot_bb(array1_matrix,array2_matrix);
      dTaudGd[ia][ib] = -1.0*dTaudGd[ia][ib]*dt_incr;
    }
  }

   //****************************************NEWTON RAPHSON CALCULATION******************************************************************************


  // Begin Newton-Raphson iterative loops.
  iNR = 0;

  converged = false;

  do {
    iNR = iNR + 1;
    // converged = true;

    NR_residual (_num_slip_sys, xs0, xm0, dt_incr, gamma_dot, F1, F_el, F_p_inv, F_p_inv_0, C,rho_SSD0,rho_SSD,rho_for_0, rho_for, bstress0, bstress, sig, tau, tau_eff, s_a, s_t, A, H, residual, sse,Prager_term,Prager_term_0,Dynamic_hardening_term,Dynamic_hardening_term_0,Static_hardening_term_0,Static_hardening_term,Energy_0,Energy,First_term_KM,First_term_KM_0,Second_term_KM,Second_term_KM_0);

    sse_ref = sse;

     //       if (sse > 0.0e0) {
     //         std::cout << "\n iNR:" << iNR;
     //         std::cout << "\n sse:" << sse << "\n";
     //       }
    
     // Begin calculation of the partial derivatives needed for the Newton-Raphson step
    
     // Calculate derivative of the mobile dislocation density, rho_for, immobile dislocation density, rho_i, backstress, bstress, w.r.t. gamma-dot-beta
       
      //  std::vector<std::vector<Real>> drhomdgb(_num_slip_sys, std::vector<Real>(_num_slip_sys));
      //  std::vector<std::vector<Real>> drhoidgb(_num_slip_sys, std::vector<Real>(_num_slip_sys));
      //  std::vector<std::vector<Real>> dbsdgb(_num_slip_sys, std::vector<Real>(_num_slip_sys));
      //  std::vector<std::vector<Real>> dldrhom(_num_slip_sys, std::vector<Real>(_num_slip_sys));
      //  std::vector<std::vector<Real>> dsadgb(_num_slip_sys, std::vector<Real>(_num_slip_sys));
      //  std::vector<std::vector<Real>> dstdgb(_num_slip_sys, std::vector<Real>(_num_slip_sys));

    std::vector<std::vector<Real>> d_rho_for_d_gamma_dot(_num_slip_sys, std::vector<Real>(_num_slip_sys));
    std::vector<std::vector<Real>> d_backstress_d_gamma_dot(_num_slip_sys, std::vector<Real>(_num_slip_sys));
    std::vector<std::vector<Real>> d_s_a_d_gamma_dot(_num_slip_sys, std::vector<Real>(_num_slip_sys));
    std::vector<std::vector<Real>> d_s_t_d_gamma_dot(_num_slip_sys, std::vector<Real>(_num_slip_sys));
    std::vector<std::vector<Real>> d_C2_d_gamma_dot(_num_slip_sys, std::vector<Real>(_num_slip_sys));
    std::vector<std::vector<Real>> d_C3_d_gamma_dot(_num_slip_sys, std::vector<Real>(_num_slip_sys));

    for(unsigned int ia = 0; ia < _num_slip_sys; ia++) {
      for(unsigned int ib = 0; ib < _num_slip_sys; ib++) {
        d_rho_for_d_gamma_dot[ia][ib] = 0.0;
        d_s_a_d_gamma_dot[ia][ib] = 0.0;
        d_backstress_d_gamma_dot[ia][ib] = 0.0;
        d_s_t_d_gamma_dot[ia][ib] = 0.0;
        d_C2_d_gamma_dot[ia][ib] = 0.0;
        d_C3_d_gamma_dot[ia][ib] = 0.0;
      }
    }

     //Calculation of constant parameters used C2 term
      Real a1=(eta_0*vol_frac*Z1)/(b_berg*lambda);
      Real a2=Z1/(b_berg*lambda);
      Real a3=Z2;
     //Forest dislocation density derivative calculation

    
     // SSD density 
    for (unsigned int ia = 0; ia < _num_slip_sys; ia++) {
      for (unsigned int ja = 0; ja < _num_slip_sys; ja++) {
        d_rho_for_d_gamma_dot[ia][ja] = (dt_incr* (k_1 * std::sqrt(rho_SSD[ia]) - k_2 * rho_SSD[ia]) * sgn(gamma_dot[ia]))/(1.0 
                   - (k_1 * std::abs(gamma_dot[ia]) * dt_incr) / (2.0* std::sqrt(rho_SSD[ia]))
                   + k_2 * std::abs(gamma_dot[ia]) * dt_incr);
        
          }
    }
     
       // Calculate derivative of the athermal slip resistance, s_a w.r.t. gamma-dot-beta.
        for(unsigned int ia = 0; ia < _num_slip_sys; ia++) {
          Real sum1 = 0.0;
          Real sum2 = 0.0;
          for(unsigned int ib = 0; ib < _num_slip_sys; ib++) {
            sum1 = 0.0; sum2 = 0.0;
            for(unsigned int ic = 0; ic < _num_slip_sys; ic++) {
                sum1 = sum1 +  A[ia][ic] * rho_SSD[ic];
                sum2 = sum2 +  A[ia][ic] * d_rho_for_d_gamma_dot[ic][ib];
            }
            d_s_a_d_gamma_dot[ia][ib] = 0.5 * k_0 * G * b_berg * sum2/sqrt(sum1);
          }
        }



        // This block computes the d(C2)/d(gamma_dot) Jacobian component
         
         // First, calculate the sum of all forest dislocation densities, as it is used repeatedly.
         Real sum_rho_for = 0.0;
         for (unsigned int ix = 0; ix < _num_slip_sys; ix++)
         {
           sum_rho_for += rho_SSD[ix];
         }
         
         // --- Safety check for the square root term ---
         if (sum_rho_for < 0.0)
         {
           sum_rho_for = 0.0; // Prevent sqrt of a negative number
         }
         Real sqrt_sum_rho = std::sqrt(sum_rho_for);
         
         // Calculate the denominator term from the formula. This term is the same for all components.
         Real denom_term = a2 + a3 * sqrt_sum_rho;
         Real denominator = denom_term * denom_term;
         
         // --- Safety check for division by zero in the denominator ---
         if (std::abs(denominator) < 1e-12)
         {
           denominator = 1e-12;
         }
         
         // --- Safety check for the term in the numerator's denominator ---
         Real num_denom_term = 2.0 * sqrt_sum_rho;
         if (std::abs(num_denom_term) < 1e-12)
         {
           num_denom_term = 1e-12;
         }
         
         // Loop over the columns of the Jacobian matrix (index beta)
         for (unsigned int ib = 0; ib < _num_slip_sys; ib++)
         {
           // Calculate the sum of the drho_for/dgb column for the current beta
           Real sum_drho = 0.0;
           for (unsigned int ix = 0; ix < _num_slip_sys; ix++)
           {
             sum_drho += d_rho_for_d_gamma_dot[ix][ib];
           }
         
           // Calculate the full numerator from the formula
           Real numerator = -a1 * a3 * (1.0 / num_denom_term) * sum_drho;
         
           // Calculate the final value for this column. It is independent of alpha (ia).
           Real value = numerator / denominator;
         
           // Assign the value to all rows in the current column
           for (unsigned int ia = 0; ia < _num_slip_sys; ia++)
           {
             d_C2_d_gamma_dot[ia][ib] = value;
           }
         }

         

           // This block computes the d(C3_prime)/d(gamma_dot) Jacobian component
           
           // Pre-calculate the sum of forest dislocation densities
           //Real sum_rho_for = 0.0;
           for (unsigned int ix = 0; ix < _num_slip_sys; ix++)
           {
             sum_rho_for += rho_SSD[ix];
           }
           
           // --- Safety check for rho_r in the denominator and exponent ---
           //Real rho_reference = rho_r; // <-- Renamed variable
           if (std::abs(rho_reference) < 1e-12)
           {
             rho_reference = 1e-12; // <-- Renamed variable
           }
           
           // Pre-calculate the constant pre-factor terms from the formula for efficiency
           Real exp_term = std::exp(-sum_rho_for / rho_reference); // <-- Renamed variable
           Real prefactor = -(r0 / rho_reference) * (1.0 - phi_s) * exp_term; // <-- Renamed variable
           
           // Loop over the columns of the Jacobian matrix (index beta)
           for (unsigned int ib = 0; ib < _num_slip_sys; ib++)
           {
             // Calculate the sum of the drho_for/dgb column for the current beta
             Real sum_drho = 0.0;
             for (unsigned int ix = 0; ix < _num_slip_sys; ix++)
             {
               sum_drho += d_rho_for_d_gamma_dot[ix][ib];
             }
           
             // Calculate the final value for this column. It is independent of alpha (ia).
             Real value = prefactor * sum_drho;
           
             // Assign the value to all rows in the current column
             for (unsigned int ia = 0; ia < _num_slip_sys; ia++)
             {
               d_C3_d_gamma_dot[ia][ib] = value;
             }
           }
           
            //Backstress derivative calculation
            // This block computes the d(Backstress)/d(gamma_dot) Jacobian component
            
         for (unsigned int ia = 0; ia < _num_slip_sys; ia++)
         {
           for (unsigned int ib = 0; ib < _num_slip_sys; ib++)
           {
             // Use a ternary operator to implement the Kronecker delta (delta_ab)
             Real delta_ab = (ia == ib ? 1.0 : 0.0);
             Real a1=(eta_0*vol_frac*Z1)/(b_berg*lambda);
             Real a2=Z1/(b_berg*lambda);
             Real a3=Z2;
             sum_rho_for += rho_SSD[ib];
             Real C2_value=a1/(a2+a3* std::sqrt(sum_rho_for));
             Real C4 = -sum_rho_for / rho_reference;
             Real C3 = r0 * (phi_s + (1.0 - phi_s) * std::exp(C4));

             // --- Calculate the Numerator ---
             Real num_term1 = C1 * delta_ab;
             Real num_term2 = -C2_value* bstress[ia] * delta_ab * sgn(gamma_dot[ia]);
             Real num_term3 = bstress[ia] * std::abs(gamma_dot[ia]) * d_C2_d_gamma_dot[ia][ib];
             Real num_term4 = bstress[ia] * d_C3_d_gamma_dot[ia][ib]; // Assuming C3 in formula is the C3' calculated before
         
             Real numerator = (num_term1 + num_term2 + num_term3 + num_term4) * dt_incr;

             // --- Calculate the Denominator ---
             Real den_term1 = delta_ab;
             Real den_term2 = C2_value * std::abs(gamma_dot[ia]) * dt_incr;
             Real den_term3 = -C3 * dt_incr; // Assuming C3 in formula is the C3' calculated before
         
             Real denominator = den_term1 + den_term2 + den_term3;
             
             
             // --- Safety check for division by zero ---
             if (std::abs(denominator) < 1e-12)
             {
               denominator = 1e-12;
             }
         
             // --- Assemble the final Jacobian component ---
              // To check the evolaution correct or not
             //d_backstress_d_gamma_dot[ia][ib] = 0;
             d_backstress_d_gamma_dot[ia][ib] = numerator / denominator;
             
           }
         }
         

         // Calculate derivative of the thermal slip resistance, s_t w.r.t. gamma-dot-beta
         for(unsigned int ia = 0; ia < _num_slip_sys; ia++) {
           for(unsigned int ib = 0; ib < _num_slip_sys; ib++) {
             d_s_t_d_gamma_dot[ia][ib] = 0.0;
           }
         }

    
    //================================Form "A-matrix" of derivatives wrt d_gamma_beta.========================================================================

    std::vector<std::vector<Real>> array3(_num_slip_sys, std::vector<Real>(_num_slip_sys));

    for(unsigned int ia = 0; ia < _num_slip_sys; ia++) {
      tau_eff[ia] = tau[ia] - bstress[ia];
      for(unsigned int ib = 0; ib < _num_slip_sys; ib++) {

       // Li Wang Jacobian matrix*************************
       //  array3[ia][ib] = -gamma_dot[ia]*(p*q*(delF0/B_k/temperature)*power(1.e0 - power((abs(tau_eff[ia])
       //      - s_a[ia])/s_t[ia], p), q - 1)*power((abs(tau_eff[ia]) - s_a[ia])/s_t[ia], p - 1))*
       //    (((dTaudGd[ia][ib] - d_backstress_d_gamma_dot[ia][ib])*sgn(tau_eff[ia]) - d_s_a_d_gamma_dot[ia][ib])/s_t[ia] - (abs(tau_eff[ia]) - s_a[ia])*d_s_t_d_gamma_dot[ia][ib]/(s_t[ia]*s_t[ia]));

       // Pin_Lu_Updated_Model*****************************
        Real r_alpha = (tau[ia] - bstress[ia]) / s_a[ia];
        Real abs_r = std::abs(r_alpha);
        Real sgn_r = sgn(tau[ia]-bstress[ia]);
        
       Real delta_ab = (ia == ib ? 1.0 : 0.0); // Kronekar delta
        array3[ia][ib] = dTaudGd[ia][ib] - d_backstress_d_gamma_dot[ia][ib] - power(std::abs(gamma_dot[ia] / gammadot0), m) * sgn(gamma_dot[ia]) * d_s_a_d_gamma_dot[ia][ib] - s_a[ia] * sgn(gamma_dot[ia]) * (m / gammadot0) * power(std::abs(gamma_dot[ia] / gammadot0), (m - 1)) * delta_ab * sgn(gamma_dot[ia]);
   
       // Changed the Jacobian Matrix *********************

       // array3[ia][ib] = dTaudGd[ia][ib]-d_backstress_d_gamma_dot[ia][ib]- d_s_a_d_gamma_dot[ia][ib]*std::pow(std::abs(gamma_dot[ia]/gammadot0),m)*(sgn(gamma_dot[ia]))- m*s_a[ia]/std::abs(gamma_dot[ia])*std::pow(std::abs(gamma_dot[ia]/gammadot0),m)*sgn(gamma_dot[ia]);
   
        }
      array3[ia][ia] = array3[ia][ia] + 1;
    }
    

    // Calculate the gradient of sse wrt gamma_dot(). Will be used later to ensure that line search is in the correct direction
    std::vector<Real> gradient(_num_slip_sys);

    for(unsigned int j = 0; j < _num_slip_sys; j++) {
      gradient[j] = 0.0;
      for(unsigned int i = 0; i < _num_slip_sys; i++) {
        gradient[j] = gradient[j] + residual[i]*array3[i][j];
      }
      gradient[j] = 2.0*gradient[j];
    }

    // Solve for increment of gamma_dot
    std::vector<Real> d_gamma_dot(_num_slip_sys);

    std::vector<std::vector<Real>> array3_inv(_num_slip_sys, std::vector<Real>(_num_slip_sys));

    MatrixTools::inverse(array3,array3_inv);

    for(unsigned int j = 0; j < _num_slip_sys; j++) {
      d_gamma_dot[j] = 0.0;
      for(unsigned int i = 0; i < _num_slip_sys; i++) {
        d_gamma_dot[j] = d_gamma_dot[j] - array3_inv[j][i]*residual[i];
      }
    }

    // Check to make sure that N-R step leads 'down hill' the sse surface
    Real sum1 = 0.0;
    for(unsigned int i = 0; i < _num_slip_sys; i++) {
      sum1 = sum1 - gradient[i]*d_gamma_dot[i];
    }

    if(sum1 > 0.0) {
      for(unsigned int i = 0; i < _num_slip_sys; i++) {
        d_gamma_dot[i] = -1*d_gamma_dot[i];
      }
    }

    // Multiply step size by two 'cause next loop will divide it by 2
    for(unsigned int k = 0; k < _num_slip_sys; k++) {
      d_gamma_dot[k] = d_gamma_dot[k]*2;
    }

    // Begin line search
    improved = false;

    N_ctr = 0;

    do {
      sse_old = sse;

      // Divide step size by 2
      for (unsigned int k = 0; k < _num_slip_sys; k++) {
        d_gamma_dot[k] = d_gamma_dot[k]/2.0;

        gamma_try[k] = gamma_dot[k] + d_gamma_dot[k];
        //gamma_try[k] = gamma_try_g[k];
      }

      NR_residual (_num_slip_sys, xs0, xm0, dt_incr, gamma_try, F1, F_el, F_p_inv, F_p_inv_0, C, rho_SSD0, rho_SSD, rho_for_0, rho_for,  bstress0, bstress, sig, tau, tau_eff, s_a, s_t, A, H, residual, sse,Prager_term,Prager_term_0,Dynamic_hardening_term,Dynamic_hardening_term_0,Static_hardening_term,Static_hardening_term_0,Energy,Energy_0,First_term_KM,First_term_KM_0,Second_term_KM,Second_term_KM_0);

      if ((sse_old <= sse_ref) && (sse >= sse_old) && (N_ctr > 0)) {
        improved = true;

        break;
      }

      N_ctr++;

    } while ((improved == false) && (N_ctr < max_loops)); // End line search

    // add d_gamma_dot to gamma_dot to get new values for this iteration
    for (unsigned int k = 0; k < _num_slip_sys; k++) {
      gamma_dot[k] = gamma_dot[k] + d_gamma_dot[k]*2.0;
    }

    // if (sse_old > tolerance) then this step has not converged
    if (sse_old <= _tol) {
      converged = true;

      // break;
    }

    // if (sse_old > sse_ref/2) then convergence is too slow and increment is divided into two subincrements
    if ((sse_old > (sse_ref/2.0)) && (converged == false)) {
      N_incr = 2*N_incr - 1;
      N_incr_total = 2*N_incr_total;

      // std::cout << "\n Total sub-increments:" << N_incr_total << "\n";
      break;
    }

  } while((iNR < 10000) && (converged == false));

  if (iNR == 10000) {
    for (unsigned int i = 0; i < 3; i++){
      for (unsigned int j = 0; j < 3; j++){
        _stress[_qp](i,j) = MooseRandom::rand()*2.0e3;
        for (unsigned int k = 0; k < 3; k++){
          for (unsigned int l = 0; l < 3; l++){
            _Jacobian_mult[_qp](i,j,k,l) = MooseRandom::rand()*2.0e3;
          }
        }
      }
    }

    throw MooseException("Crystal plasticity Newton Raphson did not converge");
    return;
  }

  // if another subincrement remains to be done, then reinitialize F_p_inv_0 and state variables.  F0 and F1 gets reinitialized back at the top of this loop
  if (N_incr < N_incr_total) {
    if (N_incr == (N_incr/2)*2) { // even
      N_incr = N_incr/2 + 1;
      N_incr_total = N_incr_total/2;
    }
    else { // odd
      N_incr = N_incr + 1;
    }

    for (unsigned int i = 0; i < 3; i++) {
      for (unsigned int j = 0; j < 3; j++) {
        F_p_inv_0(i,j) = F_p_inv(i,j);
      }
    }

    for (unsigned int k = 0; k < _num_slip_sys; k++) {
      rho_SSD0[k] = rho_SSD[k];
      bstress0[k] = bstress[k];
      Prager_term_0[k] = Prager_term[k];
      Dynamic_hardening_term_0[k] = Dynamic_hardening_term[k];
      Static_hardening_term_0[k] = Static_hardening_term[k];
      Energy_0[k]=Energy[k];
      First_term_KM_0[k]=First_term_KM[k];
      Second_term_KM_0[k]=Second_term_KM[k];
    }
  }
  else if (N_incr == N_incr_total) {
    break;
  }
  } while (N_incr_total < 10000); // End time step subincrementation

  if (N_incr_total > 10000) {
    for (unsigned int i = 0; i < 3; i++){
      for (unsigned int j = 0; j < 3; j++){
        _stress[_qp](i,j) = MooseRandom::rand()*2.0e3;
        for (unsigned int k = 0; k < 3; k++){
          for (unsigned int l = 0; l < 3; l++){
            _Jacobian_mult[_qp](i,j,k,l) = MooseRandom::rand()*2.0e3;
          }
        }
      }
    }

    throw MooseException("Crystal plasticity time step subincrement did not converge");
    return;
  }

  // Calculate the average Cauchy stress
  for (unsigned int i = 0; i < 3; i++) {
    for (unsigned int j = 0; j < 3; j++) {
      sig_avg(i,j) = sig_avg(i,j) + sig(i,j);
    }
  }

  Real F_el_matrix[3][3];
  for (unsigned int i = 0; i < 3; i++){
    for (unsigned int j = 0; j < 3; j++) {
      F_el_matrix[i][j] = F_el(i,j);
    }
  }

  // Write out Euler Angles in Bunge Convention
  aa_dot_bb(F_el_matrix,dir_cos0,dir_cos);
  bunge_angles(dir_cos,psi);

  // Calculate average elasticity tensor
  for (unsigned int i = 0; i < 3; i++) {
    for (unsigned int j = 0; j < 3; j++) {
      for (unsigned int k = 0; k < 3; k++) {
        for (unsigned int l = 0; l < 3; l++) {
          C_avg[i][j][k][l] = C[i][j][k][l];
        }
      }
    }
  }

  // Rotate xs0 and xm0 to current coordinates, called xs and xm
  F_el_inv = F_el.inverse();
  for (unsigned int n = 0; n < _num_slip_sys; n++) {
    for (unsigned int i = 0; i < 3; i++) {
      xs[i][n] = 0.0;
      xm[i][n] = 0.0;
      for (unsigned int j = 0; j < 3; j++) {
        xs[i][n] = xs[i][n] + F_el(i,j)*xs0[j][n];
        xm[i][n] = xm[i][n] + xm0[j][n]*F_el_inv(j,i);
      }
    }
  }

  // Calculate the derivative of the plastic part of the rate of deformation tensor in the current configuration wrt sigma
  for (unsigned int i = 0; i < 3; i++) {
    for (unsigned int j = 0; j < 3; j++) {
      for (unsigned int k = 0; k < 3; k++) {
        for (unsigned int l = 0; l < 3; l++) {
          ddpdsig[i][j][k][l] = 0.0e0;
        }
      }
    }
  }

 for (unsigned int i = 0; i < 3; i++) {
    for (unsigned int j = 0; j < 3; j++) {
      for (unsigned int k = 0; k < 3; k++) {
        for (unsigned int l = 0; l < 3; l++) {
          for(unsigned int n = 0; n < _num_slip_sys; n++) {
            
            tau_eff[n] = tau[n] - bstress[n];

            // --- THIS IS THE FIX ---
            // Only perform the calculation if tau_eff is not zero.
            if (std::abs(tau_eff[n]) > 1.0e-12) // Use a small tolerance
            {
              // This is your original line, now safely inside the 'if' block
              ddpdsig[i][j][k][l] = ddpdsig[i][j][k][l] + 
                  (xs[i][n]*xm[j][n] + xm[i][n]*xs[j][n]) * (((xs[k][n]*xm[l][n] + xm[k][n]*xs[l][n]) * ((1.0/m)*(gamma_dot[n]/tau_eff[n]))));
            }
            // If tau_eff[n] is zero, we add 0.0, which is correct.
            // --- END OF FIX ---
          }
        }
      }
    }
  }

  // End calculations over the integration point

  // Calculate Green strain
  E_tot = F1.transpose() * F1;
  E_tot(0,0) = E_tot(0,0) - 1.0;
  E_tot(1,1) = E_tot(1,1) - 1.0;
  E_tot(2,2) = E_tot(2,2) - 1.0;
  E_tot = 0.5 * E_tot;
 
   F_norm = F1.norm();  // returns Frobenius norm
  
  // Begin calculation of the Jacobian (the tangent stiffness matrix)

  // Store sig_avg() into sig()
  for (unsigned int i = 0; i < 3; i++) {
    for (unsigned int j = 0; j < 3; j++) {
      sig(i,j) = sig_avg(i,j);
    }
  }

  // Output stress, strain, etc. for post processing
  Real sum1;
  sum1 = 0.0;
  for(unsigned int i = 0; i < 3; i++) {
    for (unsigned int j = 0; j < 3; j++) {
      sum1 = sum1 + sig(i,j) * sig(i,j);
    }
  }

  // Calculate the inverse of F_el
  F_el_inv = F_el.inverse();

  // Scale by appropriate constants and divide by num_grains to get the average value. ALSO multiply by 'dtime' which is d(sig)/d(sig_dot)
  for (unsigned int i = 0; i < 3; i++) {
    for (unsigned int j = 0; j < 3; j++) {
      for (unsigned int k = 0; k < 3; k++) {
        for (unsigned int l = 0; l < 3; l++) {
          ddpdsig[i][j][k][l] = ddpdsig[i][j][k][l]*_dt/4.0;
        }
      }
    }
  }

  // Multiply the 4th rank elastic stiffness tensor by the derivative of the plastic part of the rate of deformation tensor wrt sig_dot
  Real array6[3][3][3][3], array4[6][6], array4_inv[6][6];
  aaaa_dot_dot_bbbb(C_avg,ddpdsig,array6);

  std::vector<std::vector<Real>> array4_matrix(6, std::vector<Real>(6));
  std::vector<std::vector<Real>> array4_inv_matrix(6, std::vector<Real>(6));

  // Add 4th rank identity tensor to array6()
  for (unsigned int i = 0; i < 3; i++) {
    for (unsigned int j = 0; j < 3; j++) {
      for (unsigned int k = 0; k < 3; k++) {
        for (unsigned int l = 0; l < 3; l++) {
          array6[i][j][k][l] = array6[i][j][k][l] + 0.5*(del[i][k]*del[j][l] + del[i][l]*del[j][k]);
        }
      }
    }
  }

  forth_to_Voigt(array6,array4);

  for (unsigned int i = 0; i < 6; i++) {
    for (unsigned int j = 0; j < 6; j++) {
      array4_matrix[i][j] = array4[i][j];
    }
  }

  // Method 2: MatrixTools::inverse
  MatrixTools::inverse(array4_matrix,array4_inv_matrix);

  for (unsigned int i = 0; i < 6; i++) {
    for (unsigned int j = 0; j < 6; j++) {
      array4_inv[i][j] = array4_inv_matrix[i][j];
    }
  }

  Voigt_to_forth(array4_inv,array6);

//   // Method 3: inverse of 4th order symmetric tensor
//   // Need to take the inverse of Array4.  Since it relates two 2nd rank tensors that are both symmetric, Array4 can be xformed to Voigt notation style in order to do the inverse, then xformed back.
//   RankFourTensor array6RFT, array6RFTinv;
//
//   for (unsigned int i = 0; i < 3; i++) {
//       for (unsigned int j = 0; j < 3; j++) {
//           for (unsigned int k = 0; k < 3; k++) {
//               for (unsigned int l = 0; l < 3; l++) {
//                   array6RFT(i,j,k,l) = array6[i][j][k][l];
//               }
//           }
//       }
//   }
//
//   array6RFTinv = array6RFT.invSymm();
//
//   for (unsigned int i = 0; i < 3; i++) {
//       for (unsigned int j = 0; j < 3; j++) {
//           for (unsigned int k = 0; k < 3; k++) {
//               for (unsigned int l = 0; l < 3; l++) {
//                   array6[i][j][k][l] = array6RFTinv(i,j,k,l);
//               }
//           }
//       }
//   }

  //  Multiply Array6 by C, the elastic stiffness matrix to finally get the Jacobian, but as a 4th rank tensor.
  aaaa_dot_dot_bbbb (array6,C_avg,ddsdde_4th);

  for (unsigned int i = 0; i < 3; i++){
    for (unsigned int j = 0; j < 3; j++) {
      for (unsigned int k = 0; k < 3; k++) {
        for (unsigned int l = 0; l < 3; l++){
          _Jacobian_mult[_qp](i,j,k,l) = ddsdde_4th[i][j][k][l];
        }
      }
    }
  }

  _stress[_qp] = sig;

  _euler_ang[_qp](0) = psi[0];
  _euler_ang[_qp](1) = psi[1];
  _euler_ang[_qp](2) = psi[2];

  // checkpoint("Storing ISVs")
  int n = -1;
  // Store dir_cos0 SDV 1-9
  for (unsigned int i = 0; i < 3; i++) {
    for (unsigned int j = 0; j < 3; j++) {
      n = n + 1;
      _state_var[_qp][n] = dir_cos0[i][j];
    }
  }

  // Store dir_cos SDV 10-18
  for (unsigned int i = 0; i < 3; i++) {
    for (unsigned int j = 0; j < 3; j++) {
      n = n + 1;
      _state_var[_qp][n] = dir_cos[i][j];
    }
  }

  // Store the elastic part of F SDV 19-27
  for (unsigned int i = 0; i < 3; i++) {
    for (unsigned int j = 0; j < 3; j++) {
      n = n + 1;
      _state_var[_qp][n] = F_el(i,j);
    }
  }

  // Store inverse of the plastic part of F SDV 28-36
  for (unsigned int i = 0; i < 3; i++) {
    for (unsigned int j = 0; j < 3; j++) {
      n = n + 1;
      _state_var[_qp][n] = F_p_inv(i,j);
    }
  }

  // Plastic strain calculations
  RankTwoTensor F_p;
  F_p = F_p_inv.inverse();

  E_p = F_p.transpose() * F_p;
  E_p(0,0) = E_p(0,0) - 1.0;
  E_p(1,1) = E_p(1,1) - 1.0;
  E_p(2,2) = E_p(2,2) - 1.0;
  E_p = 0.5 * E_p;

  // Store E_p SDV 37-45
  for (unsigned int i = 0; i<3; i++) {
    for (unsigned int j = 0; j<3; j++) {
      n = n + 1;
      _state_var[_qp][n] = E_p(i,j);
    }
  }

  // Effective strain calculations
  Real s1 = 0.0; Real s2 = 0.0;

  for (unsigned int i = 0; i<3; i++) {
    for (unsigned int j = 0; j<3; j++) {
      s1 = s1 + E_tot(i,j)*E_tot(i,j);
      s2 = s2 + E_p(i,j)*E_p(i,j);
    }
  }

  // Effective total strain
  E_eff = sqrt((2./3.) * s1);

  // Effective plastic strain
  E_p_eff = sqrt((2./3.) * s2);

  // Store E_eff SDV 46
  n = n + 1;
  _state_var[_qp][n] = E_eff;

  // Store E_p_eff SDV 47
  n = n + 1;
  _state_var[_qp][n] = E_p_eff;

  // Read backstress values SDV 48-59 
  for (unsigned int i = 0; i < _num_slip_sys; i++) {
    n = n + 1;
    _state_var[_qp][n] = bstress[i];
  }

   // backstress average calculation ISV, 60   
      for (unsigned int i = 0; i < _num_slip_sys; i++) {
      backstress_avg = backstress_avg + bstress[i];
     }
      backstress_avg = backstress_avg/_num_slip_sys;      
     // Store average backstress  SDV 60
     n = n + 1;
       _state_var[_qp][n] = backstress_avg;

   // Store GAMMA_DOT AVERAGE rates SDV-61
    
     Real gamma_dot_avg = 0.0;
     for (unsigned int i = 0; i < _num_slip_sys; i++) {
       gamma_dot_avg = gamma_dot_avg + gamma_dot[i];
     }
     gamma_dot_avg = gamma_dot_avg/_num_slip_sys;
        n = n + 1;
       _state_var[_qp][n] = gamma_dot_avg;



    //Slip resistance average value ISV 62;
     // average slip resistance value
     slip_resistance_avg=0.0;
  
     for (unsigned int i = 0; i < _num_slip_sys; i++) {
     slip_resistance_avg= slip_resistance_avg+ s_a[i];
     }
     slip_resistance_avg = slip_resistance_avg/_num_slip_sys;

     // Store average  slip resistance SDV 62
     n = n + 1;
     _state_var[_qp][n] = slip_resistance_avg;

    // Store norm of deformation gradient SDV 63
     n = n + 1;
      _state_var[_qp][n] = F_norm;

    
      // Store SSD density values SDV 64-75
      for (unsigned int i = 0; i < _num_slip_sys; i++) {
        n = n + 1;
        _state_var[_qp][n] = rho_SSD[i];
       }
      // Store average values of defect densities and loop sizes SDV 76
        rho_SSD_avg = 0.0;
        
        for (unsigned int i = 0; i < _num_slip_sys; i++) {
          rho_SSD_avg = rho_SSD_avg + rho_SSD[i];
        
        }
         rho_SSD_avg = rho_SSD_avg/_num_slip_sys;
      
      
        //SDV 76
        n = n + 1;
        _state_var[_qp][n] = rho_SSD_avg;
       

       // Read backstress first term SDV 77-88 
      for (unsigned int i = 0; i < _num_slip_sys; i++) 
      {
       n = n + 1;
       _state_var[_qp][n] = Prager_term[i];
      }

     // Read backstress second term SDV 89-100 
      for (unsigned int i = 0; i < _num_slip_sys; i++) 
      {
       n = n + 1;
       _state_var[_qp][n] = Dynamic_hardening_term[i];
      }

      // Read backstress third term SDV  101-112
      for (unsigned int i = 0; i < _num_slip_sys; i++) 
      {
       n = n + 1;
       _state_var[_qp][n] = Static_hardening_term[i];
      }

        // Store average values of backstress first term SDV 113
        Prager_term_avg = 0.0;
        
        for (unsigned int i = 0; i < _num_slip_sys; i++) {
          Prager_term_avg = Prager_term_avg + Prager_term[i];
        
        }
         Prager_term_avg = Prager_term_avg/_num_slip_sys;
      
        //SDV 113
        n = n + 1;
        _state_var[_qp][n] = Prager_term_avg;



       // Store average values of backstress second  term SDV 114
        Dynamic_hardening_term_avg= 0.0;
        
        for (unsigned int i = 0; i < _num_slip_sys; i++) {
          Dynamic_hardening_term_avg = Dynamic_hardening_term_avg + Dynamic_hardening_term[i];
        
        }
         Dynamic_hardening_term_avg = Dynamic_hardening_term_avg/_num_slip_sys;
  
        //SDV 114
        n = n + 1;
        _state_var[_qp][n] = Dynamic_hardening_term_avg;
   
       // Store average values of backstress third term  SDV 115
        Static_hardening_term_avg = 0.0;
        
        for (unsigned int i = 0; i < _num_slip_sys; i++) {
          Static_hardening_term_avg = Static_hardening_term_avg + Static_hardening_term[i];
        
        }
         Static_hardening_term_avg = Static_hardening_term_avg/_num_slip_sys;

        //SDV 115
        n = n + 1;
        _state_var[_qp][n] = Static_hardening_term_avg;


        // Read energy SDV  116-127
      for (unsigned int i = 0; i < _num_slip_sys; i++) 
      {
       n = n + 1;
       _state_var[_qp][n] = Energy[i];
      }
     
     // entropy genration rate SVD 128
      Real total_dissipation_power = 0.0;
      
      for (unsigned int i = 0; i < _num_slip_sys; i++) {
          total_dissipation_power += Energy[i];
      }

      entropy_generation_rate = total_dissipation_power / 1033.0;
      n = n + 1; 
      _state_var[_qp][n] = entropy_generation_rate;

    

       // Read  first term KM  SDV 129-140
      for (unsigned int i = 0; i < _num_slip_sys; i++) 
      {
       n = n + 1;
       _state_var[_qp][n] = First_term_KM[i];
      }   

      // Read  second term KM  SDV 129-140
      for (unsigned int i = 0; i < _num_slip_sys; i++) 
      {
       n = n + 1;
       _state_var[_qp][n] = Second_term_KM[i];
      }   


    // Store average values of backstress first term SDV 153
        First_term_KM_avg = 0.0;
        
        for (unsigned int i = 0; i < _num_slip_sys; i++) {
          First_term_KM_avg = First_term_KM_avg + First_term_KM[i];
        
        }
         First_term_KM_avg = First_term_KM_avg/_num_slip_sys;
      
        //SDV 153
        n = n + 1;
        _state_var[_qp][n] = First_term_KM_avg;



      // Store average values of KM second term SDV 154
        Second_term_KM_avg = 0.0;
        
        for (unsigned int i = 0; i < _num_slip_sys; i++) {
          Second_term_KM_avg = Second_term_KM_avg + Second_term_KM[i];
        
        }
         Second_term_KM_avg = Second_term_KM_avg/_num_slip_sys;
      
        //SDV 154
        n = n + 1;
        _state_var[_qp][n] = Second_term_KM_avg;




















  // elasticity tensor
  for (unsigned int i = 0; i < 3; ++i){
    for (unsigned int j = 0; j < 3; ++j){
      for (unsigned int k = 0; k < 3; ++k){
	    for (unsigned int l = 0; l < 3; ++l){
	      // Double check
	      _Cel_cp[_qp](i,j,k,l) = C[i][j][k][l];
	    }
      }
    }
  }

} // End computeStress()

// NR_residual()
void Pin_Lu_Full_Model::NR_residual (unsigned int num_slip_sys, std::vector<std::vector<Real>> &xs0, std::vector<std::vector<Real>> &xm0,  Real dt, std::vector<Real> gdt, RankTwoTensor F1, RankTwoTensor &F_el, RankTwoTensor &F_p_inv, RankTwoTensor F_p_inv_0, Real C[3][3][3][3],std::vector<Real> &rho_SSD0, std::vector<Real> &rho_SSD, std::vector<Real> &rho_for_0, std::vector<Real> &rho_for,  std::vector<Real> &bstress0, std::vector<Real> &bstress, RankTwoTensor &sig, std::vector<Real> &tau, std::vector<Real> &tau_eff, std::vector<Real> &s_a, std::vector<Real> &s_t, std::vector<std::vector<Real>> A, std::vector<std::vector<Real>> H, std::vector<Real> &residual, Real &sse, std::vector<Real> &Prager_term,std::vector<Real> &Prager_term_0,std::vector<Real> &Dynamic_hardening_term,std::vector<Real> &Dynamic_hardening_term_0,std::vector<Real> &Static_hardening_term,std::vector<Real> &Static_hardening_term_0,std::vector<Real> &Energy,std::vector<Real> &Energy_0,std::vector<Real> &First_term_KM,std::vector<Real> &First_term_KM_0,std::vector<Real> &Second_term_KM,std::vector<Real> &Second_term_KM_0){

  Real xL_p_inter[3][3];
  std::vector<Real> shr_g(num_slip_sys);

  // RankTwoTensor F_el;
  RankTwoTensor E_el;
  RankTwoTensor Spk2;

  std::vector<std::vector<Real>> xs(3, std::vector<Real>(num_slip_sys));
  std::vector<std::vector<Real>> xm(3, std::vector<Real>(num_slip_sys));

  // Calculate the plastic part of the velocity gradient in the intermediate configuration
  for (unsigned int i = 0; i < 3; i++) {
    for (unsigned int j = 0; j < 3; j++) {
      xL_p_inter[i][j] = 0.0;
      for (unsigned int k = 0; k < num_slip_sys; k++) {
        xL_p_inter[i][j] = xL_p_inter[i][j] + xs0[i][k]*xm0[j][k]*gdt[k];
      }
    }
  }

  // Begin calculation process of F_p_n+1 = exp(xL_p_inter*dt).F_p_n
  Real array1[3][3], array2[3][3];
  for(unsigned int i = 0; i < 3; i++) {
    for(unsigned int j = 0; j < 3; j++) {
      array1[i][j] = xL_p_inter[i][j]*dt;
    }
  }

  // Caculate Omega
  Real sum1 = 0;
  for(unsigned int i = 0; i < 3; i++) {
    for(unsigned int j = 0; j < 3; j++) {
      sum1 = sum1 + array1[i][j]*array1[i][j];
    }
  }
  Real omega = sqrt(0.5 * sum1);

  // Continue calculating intermediate stuff needed for F_p_n+1
  aa_dot_bb(array1,array1,array2);

  for(unsigned int i = 0; i < 3; i++) {
    if(omega > 0.0) { // if omega==0 then no need
    for(unsigned int j = 0; j < 3; j++) {
      array1[i][j] = array1[i][j]*sin(omega)/omega + array2[i][j]*(1.0-cos(omega))/power(omega, 2);
      }
    }
    array1[i][i] = 1.0 + array1[i][i];
  }

  // Finally multiply arrays to get F_p_inv at end of time step.
  RankTwoTensor array1_tensor, array2_tensor;

  for (unsigned int i = 0; i < 3; i++) {
    for (unsigned int j = 0; j < 3; j++) {
      array1_tensor(i,j) = array1[i][j];
    }
  }

  array2_tensor = array1_tensor.inverse();
  F_p_inv = F_p_inv_0 * array2_tensor;

  // Multiply F() by F_p_inv() to get F_el()
  F_el = F1 * F_p_inv;
  if (F_el.det() == 0.0e0) {
    F_el(0,0) = 1.0e0;
    F_el(1,1) = 1.0e0;
    F_el(2,2) = 1.0e0;
  }
  RankTwoTensor F_el_inv = F_el.inverse();

  // Rotate director vectors from intermediate configuration to the current configuration.
  for(unsigned int n = 0; n < num_slip_sys; n++) {
    for(unsigned int i = 0; i < 3; i++) {
      xs[i][n] = 0.0;
      xm[i][n] = 0.0;
      for(unsigned int j = 0; j < 3; j++) {
        xs[i][n] = xs[i][n] + F_el(i,j)*xs0[j][n];
        xm[i][n] = xm[i][n] + xm0[j][n]*F_el_inv(j,i);
      }
    }
  }

  // Calculate elastic Green Strain
  E_el = F_el.transpose() * F_el;
  E_el(0,0) = E_el(0,0) - 1.0;
  E_el(1,1) = E_el(1,1) - 1.0;
  E_el(2,2) = E_el(2,2) - 1.0;
  E_el = 0.5 * E_el;

  // Multiply the stiffness tensor by the Green strain to get the 2nd Piola Kirkhhoff stress
  Real E_el_matrix[3][3], Spk2_matrix[3][3];
  for (unsigned int i = 0; i < 3; i++) {
    for (unsigned int j = 0; j < 3; j++) {
      E_el_matrix[i][j] = E_el(i,j);
    }
  }
  aaaa_dot_dot_bb(C,E_el_matrix,Spk2_matrix);
  for (unsigned int i = 0; i < 3; i++) {
    for (unsigned int j = 0; j < 3; j++) {
      Spk2(i,j) = Spk2_matrix[i][j];
    }
  }

  // Convert from PK2 stress to Cauchy stress
  Real det = F_el.det();
  sig =  F_el * Spk2 * F_el.transpose();
  for (unsigned int i = 0; i < 3; i++) {
    for (unsigned int j = 0; j < 3; j++) {
      sig(i,j) = sig(i,j)/det;
    }
  }

  // Calculate resolved shear stress for each slip system.
  for(unsigned int k = 0; k < num_slip_sys; k++) {
    tau[k] = 0.0e0;
    for(unsigned int i = 0; i < 3; i++) {
      for(unsigned int j = 0; j < 3; j++) {
        tau[k] = tau[k] + xs[i][k]*xm[j][k]*sig(i,j);
      }
    }
  }

  
    //  calculation of ssd, backstress
       std::vector<Real> drhodt(num_slip_sys);
     //  std::vector<Real> dbstressdt(num_slip_sys);

      //calculate ssd  
      Real sum_rho = 0.0;
      for (unsigned int ib = 0; ib < num_slip_sys; ib++) {
          sum_rho = sum_rho + rho_SSD0[ib];
      }
     
      for(unsigned int ia = 0; ia < num_slip_sys; ia++) {


        First_term_KM[ia]= (k_1) * (std::sqrt(rho_SSD0[ia])) * (std::abs(gdt[ia]));
        Second_term_KM[ia]=k_2 * rho_SSD0[ia] * (std::abs(gdt[ia]));


        drhodt[ia] = (k_1) * (std::sqrt(rho_SSD0[ia])) * (std::abs(gdt[ia])) - k_2 * rho_SSD0[ia] * (std::abs(gdt[ia]));
      }
     
      for(unsigned int ia = 0; ia < num_slip_sys; ia++) {
        rho_SSD[ia] = rho_SSD0[ia] + drhodt[ia]*dt;
        if(rho_SSD[ia] < 1.0e2) {
          rho_SSD[ia] = 1.0e2;
        }
      }

        //Backstress calculation
              // --- Step 1: Pre-calculate total forest density sum ---
        Real sum_rho_for = 0.0;
        for (unsigned int ia = 0; ia < num_slip_sys; ia++)
            sum_rho_for += rho_SSD0[ia];
        
        // --- Step 2: Compute coefficients ---
        Real a1 = (eta_0 * vol_frac * Z1) / (b_berg * lambda);
        Real a2 = Z1 / (b_berg * lambda);
        Real a3 = Z2;
        
        Real C2_value = a1 / (a2 + a3 * std::sqrt(sum_rho_for));
        if (std::abs(C2_value) < 1e-12)
            C2_value = 1e-12;
        
        if (std::abs(rho_reference) < 1e-12)
            rho_reference = 1e-12;
        
        Real C4 = -sum_rho_for / rho_reference;
        Real C3 = r0 * (phi_s + (1.0 - phi_s) * std::exp(C4));
        Real C2 = C2_value;
        
        // --- Step 3: Backstress evolution for each slip system ---
        std::vector<Real> dbstressdt(num_slip_sys);
        
        for (unsigned int ia = 0; ia < num_slip_sys; ia++)
        {
            Real abs_gdot = std::abs(gdt[ia]);
            dbstressdt[ia] = C1 * gdt[ia]
                           - C2 * bstress[ia] * abs_gdot
                           + C3 * bstress[ia];
        
            // Euler Forward Update
            bstress[ia] = bstress0[ia] + dbstressdt[ia] * dt;

           //Backstress first term
            Prager_term[ia] = C1*gdt[ia];
           // second term
            Dynamic_hardening_term[ia] =  - C2 * bstress[ia] * abs_gdot;
            //Third term
            Static_hardening_term[ia] = C3 * bstress[ia];

            //energy
            Energy[ia]=(tau[ia]-bstress[ia])*gdt[ia];
            
        }

       
       // Calculate the athermal slip resistance for each slip system.
       for (unsigned int ia = 0; ia < _num_slip_sys; ia++)
       {
       // Initialize the term for the summation
       Real sum_A_rho = 0.0;
     
       // Inner loop to perform the summation over index 'xi'
       for (unsigned int ix = 0; ix < _num_slip_sys; ix++)
       {
         sum_A_rho += A[ia][ix] * rho_SSD0[ix];
       }
     
       // --- Safety check to prevent taking the square root of a negative number ---
       if (sum_A_rho < 0.0)
       {
         sum_A_rho = 0.0;
       }
     
       // Calculate the Taylor-like hardening term (the part with the square root)
       // Assuming the term 'Gb' represents the Shear Modulus * Burgers vector magnitude
       Real taylor_hardening = k_0*G * b_berg * std::sqrt(sum_A_rho);
     
       // Assemble the final slip resistance for the current slip system
       s_a[ia] = Slip_resistance_const + taylor_hardening;
      }

      // Calculate thermal slip resistance for each slip system.
      // for(unsigned int ia = 0; ia < num_slip_sys; ia++) {
      //  s_t[ia] = frictional_stress;
      //}
    
      // Calculate new slip rate for each slip system.
      for(unsigned int k = 0; k < num_slip_sys; k++) {
        Real refgdg = gammadot0;
    
        tau_eff[k] = tau[k] - bstress[k];
    
        // if(abs(tau_eff[k]) > s_a[k]) {
    
         // Flow rule for Li Wang Model
    
         // shr_g[k] = refgdg*exp((-delF0/B_k/temperature)*power(1.0 - power((abs(tau_eff[k]) - s_a[k])/s_t[k], p), q))*sgn(tau_eff[k]);
      
         // Flow rule for Pin Lu Updated Model
         
         //  Real r = std::abs(tau_eff[k]) / std::max(s_a[k], 1e-12);
          //shr_g[k] = refgdg * power(r, (1.0 / m)) * sgn(tau[k]); // see “sign” note below
          shr_g[k] = gammadot0 * power((std::abs(tau_eff[k] / s_a[k])), (1 / m)) * sgn(tau_eff[k]);
      
        //} else {
        //  shr_g[k] = 1e-18;
        // }
      }

  // Calculate residual values
  for(unsigned int k = 0; k < num_slip_sys; k++) {
    residual[k] = gdt[k] - shr_g[k];

    if (residual[k] >= 1.0e100) {
      residual[k] = 1.0e100;
    }
    if (residual[k] <= -1.0e100) {
      residual[k] = -1.0e100;
    }
  }

  // Calculate root mean square error that is to be minimized to zero
  sse = 0.0e0;
  Real gamma_dot_max = abs(gdt[0]);
  for(unsigned int k = 0; k < num_slip_sys; k++) {
    sse = sse + abs(gdt[k])*(residual[k]*residual[k]);

    if(abs(gdt[k]) > gamma_dot_max) {
      gamma_dot_max = abs(gdt[k]);
    }
  }
  if (gamma_dot_max > 1.e-18) {
    sse = sqrt(sse/gamma_dot_max)/num_slip_sys;
  }

  // std::cout << "\n sse:" << sse;
} // end NR_residual

void Pin_Lu_Full_Model::readPropsFile() {
  // MooseUtils::DelimitedFileReader reader(_propsFile);

  MooseUtils::checkFileReadable(_propsFile);
  std::ifstream file_prop;
  file_prop.open(_propsFile.c_str());

  // Assign slip system normals and slip directions for a BCC material
  for (unsigned int i = 0; i < _num_props; i++) {
    // file_prop >> _properties[_qp][i];
    if (!(file_prop >> _properties[_qp][i])) {
      mooseError("Required number of material parameters not supplied for crystal plasticity model");
    }
  }

  file_prop.close();
}

void Pin_Lu_Full_Model::assignProperties(){

  // Elastic constants
  C11 = _properties[_qp][0];           //
  C12 = _properties[_qp][1];
  C44 = _properties[_qp][2];
  G = _properties[_qp][3];
  b_berg=_properties[_qp][4];          // Burgers vector magnitude
  temperature=_properties[_qp][5];
  rho_reference=_properties[_qp][6];
  k_1=_properties[_qp][7];            
  k_2=_properties[_qp][8];
  C1=_properties[_qp][9];
  eta_0=_properties[_qp][10];
  vol_frac=_properties[_qp][11];
  Z1=_properties[_qp][12];
  Z2=_properties[_qp][13];
  lambda=_properties[_qp][14];
  delF0 = _properties[_qp][15];
  // Flow parameters
  gammadot0 = _properties[_qp][16];  // reference strain rate for Pin Lu model
  p = _properties[_qp][17];  // Shape parameter (dislocation glide)
  q = _properties[_qp][18];  // Shape parameter (dislocation glide)
  frictional_stress = _properties[_qp][19]; // Lattice frictional resistance
  Alatent = _properties[_qp][20];  // Latent hardening coefficient
  // Dislocation evolution parameters
  rho_for_zero = _properties[_qp][21];  // Initial mobile dislocation density
  B_k = _properties[_qp][22];  // Boltzmann constant
  // delF0 scaled to SI units
  //delF0 = (enthalpy_const*G*1.0e6)*(pow(b_mag,3))*1.0e-9;
  Slip_resistance_const=_properties[_qp][23];
  r0=_properties[_qp][24];
  phi_s=_properties[_qp][25];
  rho_for_avg0=_properties[_qp][26];
  m=_properties[_qp][27];
  k_0=_properties[_qp][28];
   rho_SSD_zero = _properties[_qp][29];  // Initial mobile dislocation density

}

// normalize to unit vector
void Pin_Lu_Full_Model::normalize_vector(Real *x, Real *y, Real *z) {
  Real length = sqrt(*x * *x + *y * *y + *z * *z);

  *x = *x/length;
  *y = *y/length;
  *z = *z/length;
}

// Rotate 4th order stiffness tensor
void Pin_Lu_Full_Model::rotate_4th(Real a[3][3], Real b[3][3][3][3], Real (&c)[3][3][3][3]) {
  Real d[3][3][3][3];

  for(unsigned int m = 0; m < 3; m++) {
    for(unsigned int n = 0; n < 3; n++) {
      for(unsigned int k = 0; k < 3; k++) {
        for(unsigned int l = 0; l < 3; l++) {
          d[m][n][k][l] = a[k][0]*(a[l][0]*b[m][n][0][0] + a[l][1]*b[m][n][0][1] + a[l][2]*b[m][n][0][2]) +
          a[k][1]*(a[l][0]*b[m][n][1][0] +a[l][1]*b[m][n][1][1] + a[l][2]*b[m][n][1][2]) +
          a[k][2]*(a[l][0]*b[m][n][2][0] +a[l][1]*b[m][n][2][1] + a[l][2]*b[m][n][2][2]);
        }
      }
    }
  }

  for(unsigned int i = 0; i < 3; i++) {
    for(unsigned int j = 0; j < 3; j++) {
      for(unsigned int k = 0; k < 3; k++) {
        for(unsigned int l = 0; l < 3; l++) {
          c[i][j][k][l] = a[i][0]*(a[j][0]*d[0][0][k][l] + a[j][1]*d[0][1][k][l] + a[j][2]*d[0][2][k][l]) +
          a[i][1]*(a[j][0]*d[1][0][k][l] + a[j][1]*d[1][1][k][l] + a[j][2]*d[1][2][k][l]) +
          a[i][2]*(a[j][0]*d[2][0][k][l] + a[j][1]*d[2][1][k][l] + a[j][2]*d[2][2][k][l]);
        }
      }
    }
  }
}

// 4th order tensor to Voigt notation
void Pin_Lu_Full_Model::forth_to_Voigt(Real a[3][3][3][3], Real (&b)[6][6]) {
  for(unsigned int i = 1; i <= 3; i++) {
    for(unsigned int j = i; j <= 3; j++) {
      unsigned int ia = i;
      if(i != j) {
        ia = 9 - i - j;
      }
      for(unsigned int k = 1; k <= 3; k++) {
        for(unsigned int l = k; l <= 3; l++) {
          unsigned int ib = k;
          if(k != l) {
            ib = 9 - k - l;
          }
          b[ia-1][ib-1] = a[i-1][j-1][k-1][l-1];
        }
      }
    }
  }
}

// Voigt tensor to 4th order tensor
void Pin_Lu_Full_Model::Voigt_to_forth(Real b[6][6], Real (&a)[3][3][3][3]) {
  for(unsigned int i = 1; i <= 3; i++) {
    for(unsigned int j = 1; j <= 3; j++) {
      unsigned int ia = i;
      if(i != j) {
        ia = 9 - i - j;
      }
      for(unsigned int k = 1; k <= 3; k++) {
        for(unsigned int l = 1; l <= 3; l++) {
          unsigned int ib = k;
          if(k != l) {
            ib = 9 - k - l;
          }
          a[i-1][j-1][k-1][l-1] = b[ia-1][ib-1];
          if(ia > 3) {
            a[i-1][j-1][k-1][l-1] = a[i-1][j-1][k-1][l-1]/2;
          }
          if(ib > 3) {
            a[i-1][j-1][k-1][l-1] = a[i-1][j-1][k-1][l-1]/2;
          }
        }
      }
    }
  }
}

void Pin_Lu_Full_Model::aaaa_dot_dot_bbbb(Real a[3][3][3][3], Real b[3][3][3][3], Real (&product)[3][3][3][3]) {

  for (unsigned int i = 0; i < 3; i++) {
    for (unsigned int j = 0; j < 3; j++) {
      for (unsigned int k = 0; k < 3; k++) {
        for (unsigned int l = 0; l < 3; l++) {
          product[i][j][k][l] = 0.0;
          for (unsigned int m1 = 0; m1 < 3; m1++) {
            for (unsigned int m2 = 0; m2 < 3; m2++) {
              product[i][j][k][l] = product[i][j][k][l] + a[i][j][m1][m2]*b[m1][m2][k][l];
            }
          }
        }
      }
    }
  }
}

void Pin_Lu_Full_Model::aaaa_dot_dot_bb(Real a[3][3][3][3], Real b[3][3], Real (&product)[3][3]) {

  for (unsigned int i = 0; i < 3; i++) {
    for (unsigned int j = 0; j < 3; j++) {
      product[i][j] = 0.0;
      for (unsigned int k = 0; k < 3; k++) {
        for (unsigned int l = 0; l < 3; l++) {
          product[i][j] = product[i][j] + a[i][j][k][l]*b[k][l];
        }
      }
    }
  }
}

void Pin_Lu_Full_Model::aa_dot_bb(Real a[3][3], Real b[3][3], Real (&product)[3][3]) {

  for (unsigned int i = 0; i < 3; i++) {
    for (unsigned int j = 0; j < 3; j++) {
      product[i][j] = 0.0;
      for (unsigned int k = 0; k < 3; k++) {
        product[i][j] = product[i][j] + a[i][k]*b[k][j];
      }
    }
  }
}

Real Pin_Lu_Full_Model::aa_dot_dot_bb(Real a[3][3], Real b[3][3]) {

  Real product = 0.0e0;
  for (unsigned int i = 0; i < 3; i++) {
    for (unsigned int j = 0; j < 3; j++) {
      product = product + a[i][j]*b[i][j];
    }
  }
  return product;
}

// calculate Euler angles in Bunge notation
void Pin_Lu_Full_Model::bunge_angles(Real (&array1)[3][3], Real (&psi0)[3]) {

  Real PI = 4.0*atan(1.0);
  Real phi_1, Phi, phi_2, sth;

  // From VPSC
  if (array1[2][2] >= 1.0e0) array1[2][2] = 1.0e0;
  if (array1[2][2] <= -1.0e0) array1[2][2] = -1.0e0;

  Phi = acos(array1[2][2]);

  if (abs(array1[2][2]) >= 0.99999) {
    phi_2 = 0.0;
    phi_1 = atan2(array1[0][1],array1[0][0]); // atan2(array1[0][1],array1[0][0]);
  }
  else {
    sth = sin(Phi);
    phi_2 = atan2(array1[0][2]/sth,array1[1][2]/sth);
    phi_1 = atan2(array1[2][0]/sth,-array1[2][1]/sth);
  }

  psi0[0] = 180.0e0*phi_1/PI;
  psi0[1] = 180.0e0*Phi/PI;
  psi0[2] = 180.0e0*phi_2/PI;
}

// return power (with sign)
Real Pin_Lu_Full_Model::power(Real x, Real y) {

  Real ans;

  if (x == 0.0e0) {
    if (y > 0.0e0) {
      ans = 0.0e0;
    }
    else if (y < 0.0e0) {
      ans = 1.0e300;
    }
    else {
      ans = 1.0e0;
    }
  }
  else {
    ans = y*log10(abs(x));
    if (ans > 300) {
      ans = 1.0e300;
    }
    else {
      ans = pow(10.0e0,ans);
    }
    if (x < 0.0e0) {
      ans = -ans;
    }
  }

  return ans;
}

// return signum of a real number
Real Pin_Lu_Full_Model::sgn(Real x) {

  Real ans;
  ans = 1.0;

  if (x < 0.0e0) {
    ans = -1.0;
  }

  return ans;
}

// return maximum of two scalars
Real Pin_Lu_Full_Model::max_val(Real a, Real b)
{
  if (a>b)
  {
    return a;
  }
  else
  {
    return b;
  }
}
