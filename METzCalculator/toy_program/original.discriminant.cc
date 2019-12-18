#include "TMath.h"

int main()
{
  Double_t wmass = 80.4;
    
  Double_t lmass = 0.105658367;
 
  Double_t l_x = 0.0;
  Double_t l_y = 0.0;
  Double_t l_z = 0.0;
  Double_t l_E = TMath::Sqrt(l_x*l_x + l_y*l_y + l_z*l_z + lmass*lmass);
  Double_t nu_x = 0.0;
  Double_t nu_y = 0.0;

  Double_t a = wmass*wmass - lmass*lmass + 2.0*(l_x*nu_x + l_y*nu_y);
  
  Double_t A = l_E*l_E - l_z*l_z;
  
  Double_t B = -a*l_z;
  
  Double_t C = l_E*l_E*(nu_x*nu_x + nu_y*nu_y) - 0.25*a*a;
  
  Double_t discr = B*B - 4.0*A*C;

  if (discr < 0.0)
    {
      std::cout << "bad combination.\n";

      nu_x += 1.0;
      nu_y += 1.0;

      a = wmass*wmass - lmass*lmass + 2.0*(l_x*nu_x + l_y*nu_y);
  
      A = l_E*l_E - l_z*l_z;
  
      B = -a*l_z;
  
      C = l_E*l_E*(nu_x*nu_x + nu_y*nu_y) - 0.25*a*a;
  
      discr = B*B - 4.0*A*C;
    }
