
#ifndef InverseComtponScatteringSource_h
#define InverseComtponScatteringSource_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ParticleGun.hh"
#include "globals.hh"
#include "Vector3D.hh"

#include "Randomize.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

typedef std::map<std::string, double> ConfigSection;
typedef std::map<std::string, ConfigSection> ConfigData;

struct ElectronProperties {
    // Electron's phase space properties represented with vector quantities.
    Vector3D position;  // Position vector of the electron (x, y, z).
    Vector3D beta;      // Velocity vector of the electron normalized by the speed of light (beta_x, beta_y, beta_z).
    G4double E;         // Total energy of the electron.
    G4double gamma;   // Lorentz gamma factor of the electron.

    // Default constructor
    ElectronProperties() 
        : position(Vector3D()), beta(Vector3D()), E(0.0), gamma(1.0) {}

    // Constructor to initialize the properties directly.
    ElectronProperties(const Vector3D& pos, const Vector3D& vel, G4double energy, G4double gamma)
        : position(pos), beta(vel), E(energy), gamma(gamma) {}
};

// struct ScatteredPhotonProperties {
//     // Properties of a scattered photon including its energy, position, momentum direction, and polarization represented by Stokes parameters.
//     G4double E;           // Energy of the scattered photon.
//     Vector3D position;    // Position vector where the photon is generated (x, y, z).
//     Vector3D k_vector;    // Normalized momentum vector of the photon, representing its direction.
//     Vector3D Stokes;      // Stokes parameters of the photon, representing its state of polarization.
//     Vector3D Polar_vector; // The polarization vector of the photon.

//     // Default constructor
//     ScatteredPhotonProperties()
//         : E(0.0), position(Vector3D()), k_vector(Vector3D()), Stokes(Vector3D()), Polar_vector(Vector3D()) {}

//     // Constructor to initialize the properties directly.
//     ScatteredPhotonProperties(G4double energy, const Vector3D& pos, const Vector3D& momentum, const Vector3D& stokes, const Vector3D& polar_vector)
//         : E(energy), position(pos), k_vector(momentum), Stokes(stokes), Polar_vector(polar_vector) {}
// };

struct lorentzVect {
    // 1. LorentzVect refers to the four-vector of lorentz transformation for a particle: 
    //  a. the four-momentum (energy-momentum) is (x0, x1, x2, x3) = (E/c, px, py, pz)
    //  b. the four-position (space-time) is (x0, x1, x2, x3) = (ct, x, y, z)
    // 2. In our simulation, the four-vector for a photon is set to (E, px*c, py*c, pz*c), 
    //  which is the product of (E/c, px, py, pz) with the speed of light.
    // 
    G4double x0;
    Vector3D x;
};

class G4ParticleGun;
class G4Event;

class InverseComtponScatteringSource : public G4VUserPrimaryGeneratorAction
{
    public:
        InverseComtponScatteringSource();    
        virtual ~InverseComtponScatteringSource();

        virtual void GeneratePrimaries(G4Event*);         

        const G4ParticleGun* GetParticleGun() const { return fParticleGun; }
    
    private:
        G4ParticleGun*  fParticleGun;
        G4RandFlat* RandFlat;
        G4RandGauss* RandGauss;
        G4RandExponential* RandExponential;

        // Electron beam parameters
        G4double E_e, sigE_e, nE_e;
        G4double sigX_e, sigY_e, sigZ_e, nX_e, nY_e, nZ_e;
        G4double emitX, emitY;
        G4double betaX, betaY, alphaX, alphaY, gammaX, gammaY; // twiss parameters
        G4double beamCharge;

        // Laser parameters
        G4double WL, sigWL, nWL_L, pulseEnergy;
        G4double sigX_L, sigY_L, sigZ_L;
        G4double sigX_L2, sigY_L2, sigZ_L2;
        G4double S0, S1, S2, S3, ST, tau;       
        
        // inverse Compton scattering parameters
        Vector3D pos0; // the position (X0, Y0, Z0) of the interaction point (IP)
        G4double e_i; // interaction angle
        G4int Mode, nSteps;
        G4double rangeT, dT, xsec0;

        // X/gamma
        lorentzVect ki_fV, kf_fV;

        // const value
        const G4double Pi = 3.14159265358979323846;
        const G4double twoPi = 2 * Pi;
        const G4double fourPi = 4 * Pi;
        const G4double halfPi = Pi / 2;
        const G4double Pi2 = Pi*Pi;
        const G4double reciprocal_electron_mass_c2 = 1 / electron_mass_c2;
        const G4double hc = h_Planck * c_light;
        const G4double twoPi_classic_electr_radius2 = twoPi * 2.8179403227e-15*m * 2.8179403227e-15*m;


        // Read and parse input file
        ConfigData parseConfig(const std::string& filename);

        // Determine orientation angle tau based on stokes parameter S1, S2 and S3.
        G4double Tau(const G4double, const G4double, const G4double);
        
        // Calculate the sigma of the interaction time t
        G4double sigct_cal(const G4double, const G4double, const G4double, const G4double, const G4double, const G4double);
        
        // Generate random number following Gaussian distribution with a cut tail
        G4double generateGauss(const G4double, const G4double, const G4double);

        // Sample an electron 
        ElectronProperties sampleElectron();

        // Electron drift
        void Drift(ElectronProperties&, const G4double);

        // Rotation Matrix
        Vector3D RotateY(const Vector3D&, const G4double);

        // calculate the enegy of photon based on its wavelength
        G4double wl2Ep(const G4double);

        // Calculate cross-section
        G4double xsec_ics(const G4double, const G4double, const G4double);
        G4double xsec_cal(const ElectronProperties&, Vector3D&, const G4double, const G4double);

        // Lorentz transformation
        lorentzVect lorentzTrans(const lorentzVect&, const Vector3D&, const G4double);

        // Lorentz transformation of electric field
        Vector3D lorentzTransElec(const Vector3D&, const Vector3D&, const Vector3D&, const G4double);

        // Transformation of stokes parameters
        Vector3D stokesTrans(const G4double, const G4double, const G4double, const G4double);

        // Sample the energy of the scattered photon based on the accetp-reject method
        G4double sampleE(const G4double);

        // Sample the azimethal angle psi = phi + tau of the scattered photon based on the accetp-reject method
        G4double samplePsi(const G4double a1, const G4double b1);

        // Sample the polarization vector based on the Stokes parameters of the scattered photon
        Vector3D samplePolarVec(const G4double, const G4double, const G4double, const Vector3D&, const Vector3D&);
        
        // inverse Compton scattering
        void sampleScatteredPhoton(const ElectronProperties&, const Vector3D&, G4double);
};
#endif


