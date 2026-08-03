#include "ComputeStressBase.h"
#include "EulerAngleReader.h"
#include "EBSDMeshReader.h"
#include "GrainAreaSize.h"

class Pin_Lu_Full_Model : public ComputeStressBase
{
public:
  static InputParameters validParams();

  Pin_Lu_Full_Model(const InputParameters & parameters);
  virtual ~Pin_Lu_Full_Model();

protected:
  FileName _propsFile;
  FileName _slipSysFile;

  unsigned int _num_props;
  unsigned int _num_slip_sys;
  unsigned int _num_state_vars;

  int _grainid;

  const Real _tol;
  const VariableValue & _temp;

  const EulerAngleReader * _EulerAngFileReader;
  const EBSDMeshReader * _EBSDFileReader;
  const GrainAreaSize * _GrainAreaSize;

  Real _grain_size;
  int _isEulerRadian;
  int _isEulerBunge;

  virtual void initQpStatefulProperties();
  virtual void computeQpStress();

  Real max_val(Real a,Real b);

  MaterialProperty<Point> & _euler_ang;

  // The eigenstrains
  std::vector<MaterialPropertyName> _eigenstrain_names;
  std::vector<const MaterialProperty<RankTwoTensor> *> _eigenstrains;
  std::vector<const MaterialProperty<RankTwoTensor> *> _eigenstrains_old;

  MaterialProperty<std::vector<Real> > & _state_var;
  const MaterialProperty<std::vector<Real> > & _state_var_old;
  MaterialProperty<std::vector<Real> > & _properties;
  const MaterialProperty<std::vector<Real> > & _properties_old;

  // Miller indices of slip plane normals
  MaterialProperty<std::vector<std::vector<Real> > > & _y;
  const MaterialProperty<std::vector<std::vector<Real> > > & _y_old;
  // Miller indices of slip directions
  MaterialProperty<std::vector<std::vector<Real> > > & _z;
  const MaterialProperty<std::vector<std::vector<Real> > > & _z_old;

  const MaterialProperty<RankTwoTensor> & _deformation_gradient;
  const MaterialProperty<RankTwoTensor> & _deformation_gradient_old;
  const MaterialProperty<RankTwoTensor> & _strain_increment;
  const MaterialProperty<RankTwoTensor> & _rotation_increment;
  const MaterialProperty<RankTwoTensor> & _stress_old;
  MaterialProperty<RankFourTensor> & _Cel_cp;

  // parameters/variables used for calculations
  static const int max_loops = 20;
  void readPropsFile();
  void assignProperties();
  void normalize_vector(Real*, Real*, Real*);
  void rotate_4th(Real a[3][3], Real b[3][3][3][3], Real (&c)[3][3][3][3]);
  void forth_to_Voigt(Real a[3][3][3][3], Real (&b)[6][6]);
  void Voigt_to_forth(Real b[6][6], Real (&a)[3][3][3][3]);
  void aaaa_dot_dot_bbbb(Real a[3][3][3][3], Real b[3][3][3][3], Real (&product)[3][3][3][3]);
  void aaaa_dot_dot_bb(Real a[3][3][3][3], Real b[3][3], Real (&product)[3][3]);
  void aa_dot_bb(Real a[3][3], Real b[3][3], Real (&product)[3][3]);
  Real aa_dot_dot_bb(Real a[3][3], Real b[3][3]);
  void bunge_angles(Real (&array1)[3][3], Real (&psi0)[3]);
  Real power(Real x, Real y);
  Real sgn(Real x);

  void NR_residual (unsigned int num_slip_sys, std::vector<std::vector<Real>> &xs0, std::vector<std::vector<Real>> &xm0,  Real dt, std::vector<Real> gamma_dot, RankTwoTensor F1, RankTwoTensor &F_el, RankTwoTensor &F_p_inv, RankTwoTensor F_p_inv0, Real C[3][3][3][3], std::vector<Real> &rho_SSD0, std::vector<Real> &rho_SSD,std::vector<Real> &rho_for_0, std::vector<Real> &rho_for, std::vector<Real> &bstress0, std::vector<Real> &bstress, RankTwoTensor &sig, std::vector<Real> &tau, std::vector<Real> &tau_eff, std::vector<Real> &s_a, std::vector<Real> &s_t, std::vector<std::vector<Real>> A, std::vector<std::vector<Real>> H, std::vector<Real> &residual, Real &sse, std::vector<Real> &Prager_term_0,std::vector<Real> &Prager_term,std::vector<Real> &Dynamic_hardening_term,std::vector<Real> &Dynamic_hardening_term_0,std::vector<Real> &Static_hardening_term,std::vector<Real> &Static_hardening_term_0,std::vector<Real> &Energy,std::vector<Real> &Enegry_0,std::vector<Real> &First_term_KM,std::vector<Real> &First_term_KM_0,std::vector<Real> &Second_term_KM,std::vector<Real> &Second_term_KM_0);

  Real tolerance;

  Real delF0;

  Real C11, C12, C44, G, b_berg, gammadot0, temperature,rho_reference,k_1,k_2,C1,eta_0,vol_frac,Z1,Z2,lambda,p, q,  frictional_stress, Alatent, rho_for_zero, B_k,Slip_resistance_const,r0,phi_s,rho_for_avg0,m,k_0,rho_SSD_zero;

  Real sse;
};
