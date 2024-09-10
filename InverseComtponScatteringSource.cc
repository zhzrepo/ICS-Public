#include "InverseComtponScatteringSource.hh"
#include "G4ParticleTypes.hh"
#include "G4SystemOfUnits.hh"


#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#else
#include "G4RunManager.hh"
#endif

InverseComtponScatteringSource::InverseComtponScatteringSource()
 : G4VUserPrimaryGeneratorAction(),
   fParticleGun(0)
{
    fParticleGun = new G4ParticleGun(1);

    // default particle kinematic
    fParticleGun->SetParticleDefinition(G4Gamma::Definition());

    // Try to load and apply configuration parameters
    try {
        // Initialize parameters from the configuration file
        const ConfigData config = parseConfig("input.in");

        // Lambda to safely access double values from configuration
        auto getConfigValue = [&config](const std::string& section, const std::string& key) {
            try {
                const auto& sectionMap = config.at(section);
                return sectionMap.at(key);
            } catch (const std::out_of_range&) {
                std::cerr << "Configuration for " << section << " or key " << key << " not found." << std::endl;
                std::cerr << "Please check the file \"input.in\" " << std::endl;
                throw;  // Rethrow to handle the error outside or terminate.
            }
        };
        
        // -------------------- Electron beam parameters ----------------------
        // laboratory coordinate (x, y, z)
        // (x, y, z) is a right hand coordinate system
        // z is parallel to the eletron beam propagating direction
        E_e = getConfigValue("Electron", "E_e") * MeV; // electron energy
        sigE_e = getConfigValue("Electron", "sigE_e") * E_e; // sigma of electron energy
        nE_e = getConfigValue("Electron", "nE_e"); // cutting tail of Gaussian
        
        const G4double gamma_beam = E_e * reciprocal_electron_mass_c2; // Lorentz factor, gamma
        const G4double beta_beam = sqrt(1 - 1 / (gamma_beam * gamma_beam)); // Lorentz factor, beta

        sigX_e = getConfigValue("Electron", "sigX_e") *um; // sigma size in X direction
        sigY_e = getConfigValue("Electron", "sigY_e") *um; // sigma size in Y direction
        sigZ_e = c_light * getConfigValue("Electron", "sigZ_e") *picosecond; // pulse length in Z direction
        nX_e = getConfigValue("Electron", "nX_e"); // Gaussian taill cutoff
        nY_e = getConfigValue("Electron", "nY_e"); // Gaussian taill cutoff
        nZ_e = getConfigValue("Electron", "nZ_e"); // Gaussian taill cutoff

        emitX = getConfigValue("Electron", "emitX") * mm * mrad / (beta_beam * gamma_beam); // electron beam emittance in X direction
        emitY = getConfigValue("Electron", "emitY") * mm * mrad / (beta_beam * gamma_beam); // electron beam emittance in Y direction
        betaX = sigX_e * sigX_e / emitX; // twiss parameter beta of electron in X
        betaY = sigY_e * sigY_e / emitY; // twiss parameter beta of electron in Y
        alphaX = getConfigValue("Electron", "alphaX"); // twiss parameter alpha of electron in X
        alphaY = getConfigValue("Electron", "alphaY"); // twiss parameter alpha of electron in Y
        gammaX = (1 + alphaX*alphaX) / betaX; // twiss parameter gamma of electron in X
        gammaY = (1 + alphaY*alphaY) / betaY; // twiss parameter gamma of electron in Y

        beamCharge = getConfigValue("Electron", "beamCharge") * 1e-12 * coulomb; // electron beam charge (pC)
        // ---------------------------------------------------------------------

        // -------------------- Laser parameters ----------------------
        // we assume the laser coordinate (x_L, y_L, z_L) as: 
        // 1. the laser propagates in the y-z plane of the laboratory coordinate system
        // 2. (x_L, y_L, z_L) is a right hand coordinate system
        // 3. the z_L-axis is parallel to the Laser propagating direction
        // 
        WL = getConfigValue("Laser", "WL") *nm; // wave length of laser
        sigWL = getConfigValue("Laser", "sigWL") *WL; // sigma of wave length
        nWL_L = getConfigValue("Laser", "nWL_L"); // cutting tail of Gaussian
        pulseEnergy = getConfigValue("Laser", "pulseEnergy") *joule; // sigma of wave length

        // Z0 axis : parallel to the Laser direction
        // X0-Y0 : perpendicular to the Z_L axis
        sigX_L = getConfigValue("Laser", "sigX_L") *um; // sigma size of the waist in X_L direction
        sigY_L = getConfigValue("Laser", "sigY_L") *um; // sigma size of the waist in Y_L direction 
        sigZ_L = c_light * getConfigValue("Laser", "sigZ_L") *picosecond; // pulse length in Z0 direction

        S0 = 1; // stokes parameter S0, which is usually set to 1 unless specific settings are required
        S1 = getConfigValue("Laser", "S1"); // stokes parameter S1
        S2 = getConfigValue("Laser", "S2"); // stokes parameter S2
        S3 = getConfigValue("Laser", "S3"); // stokes parameter S3
        // -------------------------------------------------------------

        // --------------------- inverse Compton scattering ----------------------
        pos0 = Vector3D(getConfigValue("Setup", "X0"), getConfigValue("Setup", "Y0"), getConfigValue("Setup", "Z0")); // the position of the IP
        e_i = getConfigValue("Setup", "e_i") * deg; // colliding angle between the laser and the electron beam
        Mode = getConfigValue("Setup", "Mode"); // simulating mode
        nSteps = getConfigValue("Setup", "nSteps"); // Since the time of interaction is divided into equal steps, nSteps is the number of the steps
        // -----------------------------------------------------------------------

        // ---------------------- Pre calculation  -----------------------
        // Check the stokes parameters
        if ((S1*S1 + S2*S2 + S3*S3) > 1) {
            G4cout << "(S1*S1 + S2*S2 + S3*S3) > 1" << G4endl;
            G4cout << "Please check the file \"input.in\". " << G4endl;
            exit(EXIT_FAILURE);
        }

        // Check the colliding angle
        if (e_i < 5*deg) {
            G4cout << "e_i < 5*deg" << G4endl;
            G4cout << "Please check the file \"input.in\". " << G4endl;
            exit(EXIT_FAILURE);
        }

        // Set sigxl^2, sigyl^2, sigzl^2, ST and tau
        sigX_L2 = sigX_L * sigX_L;
        sigY_L2 = sigY_L * sigY_L;
        sigZ_L2 = sigZ_L * sigZ_L;
        ST = std::sqrt(S1 * S1 + S3 * S3);
        tau = Tau(S1, S2, S3); // Determine the orientation angle tau based on the degree of circular polarization

        // Calculate the range of time t
        rangeT = 5 * sigct_cal(sigX_e, sigZ_e, sigX_L, sigZ_L, beta_beam, e_i);
        nSteps += nSteps % 2 == 0;
        dT = 2 * rangeT / (nSteps - 1);

        // Determine the simulating mode and set relevant parameters.
        if (Mode == 1) {
            const G4double rayX = fourPi * sigX_L2 / (WL + nWL_L * sigWL);
            const G4double rayY = fourPi * sigY_L2 / (WL + nWL_L * sigWL);
            const G4double cx2 = nX_e * nX_e * sigX_e * sigX_e * nZ_e * nZ_e * sigZ_e * sigZ_e / (rayX * rayX + nZ_e * nZ_e * sigZ_e * sigZ_e);
            const G4double cy2 = nY_e * nY_e * sigY_e * sigY_e * nZ_e * nZ_e * sigZ_e * sigZ_e / (rayY * rayY + nZ_e * nZ_e * sigZ_e * sigZ_e);
            const G4double alpha1 = std::atan(std::sqrt(gammaX * emitX * nX_e * nX_e + gammaY * emitY * nY_e * nY_e));
            const G4double alpha2 = std::atan(std::sqrt(cx2 + cy2));
            const G4double e_imax = (e_i + alpha1 + alpha2 > 180 * deg) ? Pi : e_i + alpha1 + alpha2;
            const G4double e_imin = (e_i - alpha1 - alpha2 > 0 * deg) ? e_i - alpha1 - alpha2 : 0;
            const G4double gamma_max = (E_e + nE_e * sigE_e) * reciprocal_electron_mass_c2;
            const G4double beta_max = std::sqrt(1 - 1 / gamma_max / gamma_max);
            xsec0 = sigX_L * sigY_L / (1 - std::cos(e_imax)) / xsec_ics(WL - nWL_L * sigWL, 1 - beta_max * std::cos(e_imin), gamma_max);
        } else if (Mode == 2) {
            const G4double ep = wl2Ep(WL);
            xsec0 = pulseEnergy / ep / sigZ_L  * (0.25 / Pi / Pi) * dT * std::sqrt(twoPi);

            if (xsec0 * 2 * 0.665 * barn > 1) {
                G4cout << "The probability is bigger than 1. Please increase the number of nSteps." << G4endl;
            }

        } else {
            std::cerr << "Please set the simulating Mode." << std::endl;
            exit(EXIT_FAILURE);
        }
        // ---------------------------------------------------------------
    } catch (const std::out_of_range& e) {
        std::cerr << "Configuration access error: " << e.what() << std::endl;
        exit(EXIT_FAILURE); // Exit if configuration is critical
    } catch (const std::exception& e) {
        std::cerr << "Unexpected error during initialization: " << e.what() << std::endl;
        exit(EXIT_FAILURE); // Exit if other unexpected errors occur
    }
}

InverseComtponScatteringSource::~InverseComtponScatteringSource()
{
    // This block executes only if the thread ID is 0 (the master thread in a multithreaded context) and Mode is set to 2.
    if (G4Threading::G4GetThreadId() == 0 && Mode == 2) {
        G4RunManager* runManager = G4RunManager::GetRunManager();
        // Outputs the equivalent number of real particles each macro-particle represents.
        // The calculation uses the total beam charge and the number of events to be processed.
        G4cout << "One macro-particle represents " << beamCharge / runManager->GetNumberOfEventsToBeProcessed() << " real particles." << G4endl;
    }
    delete fParticleGun;
}

ConfigData InverseComtponScatteringSource::parseConfig(const std::string& filename)
{
    ConfigData config;
    std::ifstream file(filename);
    std::string line, currentSection;

    if (!file.is_open()) {
        std::cerr << "Unable to open file: " << filename << std::endl;
        exit(1);
    }
    
    while (getline(file, line)) {
        // Remove comments
        size_t commentPos = line.find("//");
        if (commentPos != std::string::npos) {
            line = line.substr(0, commentPos);
        }

        // Trim whitespace from the beginning and end
        line.erase(0, line.find_first_not_of(" \t\n\r\f\v"));
        line.erase(line.find_last_not_of(" \t\n\r\f\v") + 1);

        // Check if the line is empty after trimming and removing comments
        if (line.empty()) {
            continue;
        }

        // Check for and handle section headers
        if (line.front() == '[' && line.back() == ']') {
            currentSection = line.substr(1, line.size() - 2);
            continue;
        }

        // Parse key-value pairs
        std::istringstream iss(line);
        std::string key, value;
        if (getline(iss, key, '=') && getline(iss, value)) {
            key.erase(key.find_last_not_of(" \t\n\r\f\v") + 1);
            value.erase(0, value.find_first_not_of(" \t\n\r\f\v"));
            if (!key.empty() && !value.empty()) {
                try {
                    config[currentSection][key] = std::stod(value);
                } catch (std::exception& e) {
                    std::cerr << "Failed to convert " << value << " to double for key " << key << " in section " << currentSection << std::endl;
                }
            }
        }
    }

    file.close();
    return config;
}

G4double InverseComtponScatteringSource::Tau(const G4double s1, const G4double s2, const G4double s3)
{
    const G4double sT = std::sqrt(s1 * s1 + s3 * s3);
    if (!(std::abs(s2) == 1)) {
        // If the polarization is not circular, calculate tau as half the angle from the arctangent of normalized S1 and S3
        // Normalize the angle to [0, π)
        return std::fmod(atan2(s1, s3) + twopi, twopi) / 2;
    } else {
        // If the polarization is circular, set tau to 0
        return 0;
    }
}

G4double InverseComtponScatteringSource::sigct_cal(const G4double sigxe, const G4double sigze, const G4double sigxl, const G4double sigzl, const G4double beta, const G4double alpha)
{
    const G4double cosAlpha = std::cos(alpha);
    const G4double sinAlpha = std::sin(alpha);

    // Compute square of trigonometric functions
    const G4double cosAlphaSquared = cosAlpha * cosAlpha;
    const G4double sinAlphaSquared = sinAlpha * sinAlpha;

    // Calculating terms used in the numerator
    const double sigxeSquared = sigxe * sigxe;
    const double sigzeSquared = sigze * sigze;
    const double sigxlSquared = sigxl * sigxl;
    const double sigzlSquared = sigzl * sigzl;

    const double term1 = sigxlSquared * (1 - beta * cosAlpha) * (1 - beta * cosAlpha) + sigxeSquared * (beta - cosAlpha) * (beta - cosAlpha);
    const double term2 = (beta * beta * sigzlSquared + sigzeSquared) * sinAlphaSquared;
    const double denominator = term1 + term2;
    // Check for zero denominator
    if (denominator == 0) {
        throw std::runtime_error("Division by zero in sigct_cal.");
    }

    const double term3 =  (1 - beta) * (1 - beta) * (sigzlSquared * sigxlSquared + sigzlSquared * sigxeSquared * cosAlphaSquared + sigxeSquared * sigxlSquared * sinAlphaSquared);
    const double term4 = sigzeSquared * (1 - cosAlpha) * (1 - cosAlpha) * (sigxlSquared + sigxeSquared) + sigzlSquared * sigzeSquared * sinAlphaSquared;
    const double numerator = term3 + term4;

    // Compute the final result
    return std::sqrt(numerator / denominator);
}

G4double InverseComtponScatteringSource::generateGauss(const G4double u, const G4double sig, const G4double n)
{
    // Function to generate a Gaussian-distributed random value.
    //
    // Parameters:
    //   u   - Mean of the Gaussian distribution.
    //   sig - Standard deviation of the Gaussian distribution.
    //   n   - Cutoff multiplier to limit the range of the distribution to within ±n standard deviations.

    G4double r;
    do {
        r = RandGauss->shoot();
    } while (std::abs(r) > n);
    return u + sig * r;
}

ElectronProperties InverseComtponScatteringSource::sampleElectron()
{
    // Function to sample an electron from the phase space of the electrons.
    // The return type of ElectronProperties can be found in InverseComtponScatteringSource.hh

    // const G4double twoPi = 2 * CLHEP::pi;

    // Lambda for generating phase-space coordinates [S, p_S]
    auto generatePhaseSpaceCoordinates = [this](double emit, double beta, double alpha, double n, auto& randExp, auto& randFlat) {
        G4double u, phi, S, p_S;
        do {
            u = randExp->shoot(1.0);
        } while (u > n);
        phi = randFlat->shoot(0., twoPi);
        S = sqrt(2 * u * emit * beta) * cos(phi);
        p_S = -sqrt(2 * u * emit / beta) * (alpha * cos(phi) + sin(phi));
        return std::make_pair(S, p_S);
    };

    // Generating x, xp, y, yp, z, and E
    auto [x, xp] = generatePhaseSpaceCoordinates(emitX, betaX, alphaX, nX_e * nX_e / 2, RandExponential, RandFlat);
    auto [y, yp] = generatePhaseSpaceCoordinates(emitY, betaY, alphaY, nY_e * nY_e / 2, RandExponential, RandFlat);
    const G4double z = generateGauss(0, sigZ_e, nZ_e);
    const G4double E = generateGauss(E_e, sigE_e, nE_e);

    // Normalize the velocity using the speed of light
    const G4double gamma = E * reciprocal_electron_mass_c2;
    const G4double beta = sqrt(1 - 1 / (gamma * gamma));
    const G4double norm = 1 / sqrt(xp * xp + yp * yp + 1);
    
    return ElectronProperties(Vector3D(x, y, z), Vector3D(beta * xp * norm, beta * yp * norm, beta * norm), E, gamma);
}

void InverseComtponScatteringSource::Drift(ElectronProperties& electron, const G4double dt)
{
    electron.position = electron.position + dt * electron.beta;
}

Vector3D InverseComtponScatteringSource::RotateY(const Vector3D& vect, const G4double phi)
{
    // Performs a counterclockwise rotation of the given vector 'vect' around the Y-axis by the angle 'phi'.
    // Frames and rotation:
    // - All frames are the right-handed Cartesian coordinate system.
    // - The rotation of the vector is counterclockwise when viewed from the top of the Y-axis.
    // - Note: The counterclockwise rotation of the vector corresponds to a clockwise rotation of the coordinate system.
    //   Thus, the angle 'phi' is used as is for matrix application, assuming counterclockwise vector rotation.

    // Parameters:
    // vect: The vector to be rotated, passed by constant reference to avoid unnecessary copying.
    // phi: The rotation angle.

    // Output:
    // Returns the Vector3D containing the new coordinates of the vector after rotation around the Y-axis.
    // The x and z coordinates are adjusted based on the rotation angle, while the y coordinate remains unchanged.

    const G4double sinphi = std::sin(phi);  // Sine of the rotation angle
    const G4double cosphi = std::cos(phi);  // Cosine of the rotation angle

    // Compute new coordinates by applying the rotation matrix:
    return Vector3D(
        cosphi * vect.x + sinphi * vect.z,  // New x-coordinate after rotation
        vect.y,                             // y-coordinate remains unchanged
        -sinphi * vect.x + cosphi * vect.z  // New z-coordinate after rotation
    );
}

G4double InverseComtponScatteringSource::wl2Ep(const G4double wl)
{
    return hc / wl;
}

G4double InverseComtponScatteringSource::xsec_ics(const G4double wl, const G4double one_minus_beta_dot_x, const G4double gamma) 
{
    const G4double X = 2 * gamma * wl2Ep(wl) * one_minus_beta_dot_x * reciprocal_electron_mass_c2;
    const G4double reciprocal_X = 1 / X;
    const G4double reciprocal_Xp1 = 1 - 1 / ((X + 1) * (X + 1));
    if (X > 0.01) {
        return twoPi_classic_electr_radius2 * reciprocal_X * ((1 - 4 * reciprocal_X - 8 * reciprocal_X * reciprocal_X) * std::log(1 + X) + 8 * reciprocal_X + 0.5 * reciprocal_Xp1);
    } else {
        return twoPi_classic_electr_radius2 * 1.33333333333333 * (1 - X);
    }
}

G4double InverseComtponScatteringSource::xsec_cal(const ElectronProperties& electron, Vector3D& k_vector, const G4double wl, const G4double tt) 
{
    // Transport the electron position (x, y, z) from the lab frame to the laser frame (x_L, y_L, z_L).
    // We define this rotation, which is from the lab frame to laser frame, is counterclockwise.
    // Based on the defination of RotateY, the input variable here is -e_i.
    Vector3D e_pos = electron.position;
    auto [x_L, y_L, z_L] = RotateY(e_pos, -e_i);
    
    // Calculate the Rayleigh length for the laser in the x and y directions
    const G4double rayl_X = fourPi * sigX_L2 / wl;
    const G4double rayl_Y = fourPi * sigY_L2 / wl;

    // Adjust the sigma values based on the position in the laser coordinate
    const G4double sigX2 = sigX_L2 * (1 + z_L * z_L / (rayl_X * rayl_X));
    const G4double sigY2 = sigY_L2 * (1 + z_L * z_L / (rayl_Y * rayl_Y));
    const G4double sigX = std::sqrt(sigX2);
    const G4double sigY = std::sqrt(sigY2);

    // Calculate the normalized k-vector components in the laser frame
    const G4double c1 = x_L * z_L / (rayl_X * rayl_X + z_L * z_L);
    const G4double c2 = y_L * z_L / (rayl_Y * rayl_Y + z_L * z_L);
    const G4double cFactor = 1 / std::sqrt(1 + c1 * c1 + c2 * c2);
    const G4double kx_L = c1 * cFactor;
    const G4double ky_L = c2 * cFactor;
    const G4double kz_L = cFactor;

    // Rotate k-vector to the lab frame
    k_vector = RotateY(Vector3D(kx_L, ky_L, kz_L), e_i); // laser-frame to lab-frame, 

    // Compute the exponential factors for the normalized probability density function
    const G4double u1 = x_L * x_L / (2 * sigX2) + y_L * y_L / (2 * sigY2);
    const G4double u2 = (z_L - tt) * (z_L - tt) / (2 * sigZ_L2);
    const G4double xsec1 = 1 - k_vector.dot(electron.beta); // 1 - k-vector * beta-vector
    const G4double xsec2 = exp(-u1 - u2) / sigX / sigY;
    return xsec0 * xsec1 * xsec2 * xsec_ics(wl, xsec1, electron.gamma);
}

lorentzVect InverseComtponScatteringSource::lorentzTrans(const lorentzVect& X, const Vector3D& beta, const G4double gamma)
{
    // Performs a Lorentz transformation on a four-vector X given the velocity vector beta.
    // The transformation accounts for relativistic effects when changing reference frames.

    // Parameters:
    // X: The four-vector to be transformed, represented as (time component, x, y, z) or equivalently as (E, px*c, py*c, pz*c),
    //    where E is energy, px, py, pz are the momentum components, and c is the speed of light.
    // beta: The velocity vector (vx, vy, vz) of the new frame relative to the old frame, normalized by the speed of light.
    // gamma: The Lorentz factor, calculated from the magnitude of beta (gamma = 1.0 / sqrt(1.0 - beta^2)).

 
    // Output:
    // Returns the Lorentz-transformed four-vector.

    // Formula Source:
    // The transformation is based on the Lorentz transformation matrix, B(v), which is derived from the theory of Special Relativity.
    // x0' = gamma(x0 - beta .* x); 
    // x' = x + ((gamma - 1) * (beta .* x) / beta^2 - gamma * x0) * beta 
    //    = x + (gamma^2 / (gamma + 1) * (beta .* x) - gamma * x0) * beta
    // For more detailed formulas and theory, see: https://en.wikipedia.org/wiki/Lorentz_transformation#Matrix_forms

    // Extract the components of the four-vector and the velocity vector.
    auto [x0, x] = X;

    // Calculate the dot product of beta and the spatial(momentum) components of X.
    const G4double beta_dot_X = beta.dot(x);

    // Perform the Lorentz transformation for the time component.
    const G4double x0_prime = gamma * (x0 - beta_dot_X);

    // Perform the Lorentz transformation for the spatial components.
    const G4double factor = gamma * gamma / (gamma + 1) * beta_dot_X - gamma * x0;
    const Vector3D x_prime = x + factor * beta;

    // Return the transformed four-vector.
    return lorentzVect{x0_prime, x_prime};
}

Vector3D InverseComtponScatteringSource::lorentzTransElec(const Vector3D& E_vec, const Vector3D& B_vec, const Vector3D& beta, const G4double gamma) 
{
    const Vector3D unitBeta = beta.normalize();
    return (E_vec + beta.cross(B_vec)) * gamma - unitBeta * (gamma - 1) * (E_vec.dot(unitBeta));
}

Vector3D InverseComtponScatteringSource::stokesTrans(const G4double s1, const G4double s2, const G4double s3, const G4double phi)
{
    // Performs the transformation of the Stokes parameters when the coordinate
    // system of polarization is rotated.
    
    // Parameters:
    // s1, s2, s3: The Stokes parameters before the rotation.
    // phi: The angle of the counterclockwise rotation.

    // Output:
    // Returns the Stokes parameters after the rotation. 

    const G4double cos_2phi = std::cos(2 * phi);
    const G4double sin_2phi = std::sin(2 * phi);

    return Vector3D(s1 * cos_2phi - s3 * sin_2phi,
                    s2,
                    s1 * sin_2phi + s3 * cos_2phi);
}


G4double InverseComtponScatteringSource::sampleE(const G4double Ei)
{
    // Function to sample energy based on rejection sampling.
    // Parameters:
    //   Ei - The initial energy for which the distribution parameters are to be calculated.
    
    // Calculate constants a, f1 (a^2), f2 (a^2 + 2*a), and f3 (2a^2 + 2a - 1) for the equations.
    const G4double a = electron_mass_c2 / Ei; // Ratio of electron mass energy to photon energy, defining the energy scale.
    const G4double f1 = a * a; 
    const G4double f2 = f1 + 2 * a; 
    const G4double f3 = f1 + f2 - 1;

    // Define the Probability density function (PDF) of the energy for rejection sampling.
    // PDF: f(t) = a^2/t^2 + t + a^2 + 2a - (2a^2 + 2a - 1)/t
    auto f = [f1, f2, f3](double t) { return f1 / (t * t) + t + f2 - f3 / t; };

    // Calculate the minimum 't' and the range over which 't' should be sampled.
    const G4double t_min = a / (a + 2);
    const G4double t_range = 2 / (a + 2);

    // Determine the maximum value of the function to scale the acceptance probability.
    const G4double f_max = 2 + 4 / f2;

    G4double t; // Variable to store the sampled value.

    // Perform rejection sampling to find a value 't' that fits the defined distribution.
    do {
        t = RandFlat->shoot(0.0, 1.0) * t_range + t_min; // Randomly sample a value 't' within the range from 't_min' to 't_min' + 't_range'.
    } while ((f_max * RandFlat->shoot(0.0, 1.0)) > f(t)); // Check if the sampled 't' is accepted based on the ratio of its function value to the maximum.

    // Scale the accepted 't' value by the initial energy to get the sampled energy.
    return t * Ei;
}

G4double InverseComtponScatteringSource::samplePsi(const G4double a1, const G4double b1)
{
    // Function to sample an angle psi based on a modified cosine function using rejection sampling.
    // Parameters:
    //   a1 - Adjusts the baseline height of the cosine function.
    //   b1 - Modifies the amplitude of the cosine function.

    // Define the function h(p) to be used for rejection sampling,
    // h(p) = a1 - b1 * cos(2p)
    auto h = [a1, b1](double p) { return a1 - b1 * std::cos(2*p); };

    // Set the minimum and maximum values for angle p, which spans the full circle from 0 to 2*PI.
    const G4double p_min = 0;
    const G4double p_range = twoPi;

    // Determine the maximum value of function h, used to scale the acceptance probability.
    const G4double h_max = a1 + b1;

    G4double p; // Variable to store the sampled angle.

    // Perform rejection sampling to find an angle 'p' that fits the defined function h(p).
    do {
        // Generate a random 't' value within the specified range
        p = RandFlat->shoot(0.0, 1.0) * p_range; // Randomly generate an angle 'p' within the full circular range [0, 2*PI).
    } while ((h_max * RandFlat->shoot(0.0, 1.0)) > h(p)); // Check if the sampled angle 'p' is accepted based on the ratio of its function value to the maximum.

    // Return the accepted angle, which correctly follows the distribution defined by h.
    return p;
}

Vector3D InverseComtponScatteringSource::samplePolarVec(const G4double s1, const G4double s2, const G4double s3, const Vector3D& ex, const Vector3D& ey)
{   
    // Calculate s0 of the Stokes parameters
    G4double s0 = sqrt(s1 * s1 + s2 * s2 + s3 * s3);

    // Compute the polarization angle tau_polar
    G4double tau_polar = Tau(s1, s2, s3);

    // Compute sin_2chi = sin(2*chi), where chi is the polarization azimuthal angle
    G4double sin_2chi = s2 / s0;
    G4double cos_4chi = 1 - 2 * sin_2chi * sin_2chi;

    // Determine the correct form of phi based on the value of cos(4*chi)
    G4double phi;
    if (cos_4chi > 0) {
        phi = samplePsi(1, cos_4chi) + tau_polar + halfPi;
    } else {
        phi = samplePsi(1, -cos_4chi) + tau_polar;
    }

    // Return the polarization vector as a linear combination of the base vectors ex and ey
    return std::cos(phi) * ex + std::sin(phi) * ey;    
}

void InverseComtponScatteringSource::sampleScatteredPhoton(const ElectronProperties& electron, const Vector3D& ki, G4double wl)
{
    // ------------------------------------- Sample scattered photon -------------------------------------
    // 1. Determine the four-vector and the electric vector of the incident photon
    // ki_fV (E, px*c, py*c, pz*c) is the four-vector of incident photon, which is the product
    // of (E/c, px, py, pz) with the speed of light. "fV" stands for "four-vector".
    const G4double E_photon = wl2Ep(wl);
    ki_fV = lorentzVect{E_photon, ki * E_photon};
    const Vector3D ex = Vector3D(ki.z, 0, -ki.x).normalize(); // the parallel-component of the electric field

    // 2. Transport the incident photon from the lab frame to the electron rest frame
    lorentzVect ki_fV_prime = lorentzTrans(ki_fV, electron.beta, electron.gamma);
    Vector3D ExPrime_vec = lorentzTransElec(ex, ki.cross(ex), electron.beta, electron.gamma);

    // 3: Sample the energy and the azimethal angle phi of the scattered photon, 
    //    determine the four-vector of the scattered photon in the electron rest frame,
    //    and determine the Stokes parameters of the scattered photon in the Fano reference system

    // Sample the energy Ef of the scattered photon, and calculate the scattered angle theta
    const G4double Ef = sampleE(ki_fV_prime.x0); 
    const G4double cos_e = 1 - (1 / Ef - 1 / ki_fV_prime.x0) * electron_mass_c2;
    const G4double sin_e2 = 1 - cos_e * cos_e;
    const G4double sin_e = std::sqrt(sin_e2);

    // Sample the azimethal angle phi, and calculate the Stokes parameters of the scattered photon
    const G4double tFactor = Ef / ki_fV_prime.x0 + ki_fV_prime.x0 / Ef;
    G4double psi, phi, F, PF1, PF2, PF3;
    if (ST != 0) {
        psi = samplePsi(tFactor - sin_e2, ST * sin_e2);
        phi = psi + tau;

        F = tFactor - (1 + ST * std::cos(2 * psi)) * sin_e2;
        PF1 = 2 * cos_e * ST * std::sin(2 * psi) / F;
        PF2 = tFactor * cos_e * S2 / F;
        PF3 = (sin_e2 - (2 - sin_e2) * ST * std::cos(2 * psi)) / F;
    } else {
        // Circular polarization
        phi = RandFlat->shoot(0.0, 1.0) * twoPi;

        F = tFactor - sin_e2;
        PF1 = 0;
        PF2 = tFactor * cos_e / F * S2;
        PF3 = sin_e2 / F;
    }

    // Based on the theta and phi, determine the momentum direction of the scattered photon.
    const Vector3D ez_prime = (ki_fV_prime.x).normalize(); // ez': the direction of the incident photon
    const Vector3D ex_prime = ExPrime_vec.normalize(); // ex': the direction of the electric component
    const Vector3D ey_prime = ez_prime.cross(ex_prime); // ey': cross product of Z'_L and X'_L, the same as the lorentz transformation of the Ey_vec.
    const Vector3D kf_prime = ez_prime * cos_e + ex_prime * sin_e * std::cos(phi) + ey_prime * sin_e * std::sin(phi); // Ensure this is a unit vector.

    // 4. Execute the following transformations:
    //    - Transform the four-vector of the scattered photon from the electron rest frame to the laboratory frame.
    //    - Convert the Stokes parameters of the scattered photon from the Fano reference system to the coordinate 
    //      system in the laboratory frame that describes the polarization of the scattered photon.
    kf_fV = lorentzTrans(lorentzVect{Ef, kf_prime * Ef}, electron.beta * (-1), electron.gamma);

    const Vector3D ex_Fano_prime = (ez_prime.cross(kf_prime)).normalize();
    const Vector3D ey_Fano_prime = kf_prime.cross(ex_Fano_prime);
    const Vector3D ex_Fano = lorentzTransElec(ex_Fano_prime, ey_Fano_prime, electron.beta * (-1), electron.gamma);
    const Vector3D kf = (kf_fV.x).normalize();
    const Vector3D ex_Scat = (Vector3D(kf.z, 0, -kf.x)).normalize();
    const Vector3D ey_Scat = kf.cross(ex_Scat);
    // Alternative polarization coordinate system:
    // 1. 
    // const Vector3D ex_lab = ex;
    // const Vector3D ex_Scat = (ex_lab - ex_lab.dot(kf) * kf).normalize();
    // const Vector3D ey_Scat = kf.cross(ex_Scat);
    // 2. 
    // const Vector3D ex_Scat = Vector3D(1, 0, 0);
    // const Vector3D ey_Scat = Vector3D(0, 1, 0);
    
    const auto& [stokes1, stokes2, stokes3] = stokesTrans(PF1, PF2, PF3, -atan2(ex_Fano.dot(ey_Scat), ex_Fano.dot(ex_Scat)));
    // ---------------------------------------------------------------------------------------------------
    
    const auto& [polarX, polarY, polarZ] = samplePolarVec(stokes1, stokes2, stokes3, ex_Scat, ey_Scat);
    fParticleGun->SetParticleEnergy(kf_fV.x0);
    fParticleGun->SetParticleMomentumDirection(G4ThreeVector(kf.x, kf.y, kf.z));
    fParticleGun->SetParticlePolarization(G4ThreeVector(stokes1, stokes2, stokes3));
    // fParticleGun->SetParticlePolarization(G4ThreeVector(polarX, polarY, polarZ));
}

void InverseComtponScatteringSource::GeneratePrimaries(G4Event* anEvent)
{
    ElectronProperties electron;
    Vector3D k_vector;
    G4double wl; // wave length of the photon
    G4double tt; // the time of the interaction
    G4double r, xsec;
    if (Mode == 1) {
        do {
            // Sample an electron, the wave length of the photon and the interaction time.
            electron = sampleElectron();
            wl = generateGauss(WL, sigWL, nWL_L);
            // tt = generateGauss(0.0, sig_ct, n_ct); // unit: mm
            tt = RandFlat->shoot(-rangeT, rangeT);
            Drift(electron, tt);
            r = RandFlat->shoot(0.0, 1.0);
            xsec = xsec_cal(electron, k_vector, wl, tt);
        } while (r > xsec);
        sampleScatteredPhoton(electron, k_vector, wl);
        const auto& [posX, posY, posZ] = electron.position + pos0;
        fParticleGun->SetParticlePosition(G4ThreeVector(posX, posY, posZ));
        fParticleGun->GeneratePrimaryVertex(anEvent);
    } else if (Mode == 2) {
        electron = sampleElectron();
        wl = generateGauss(WL, sigWL, nWL_L);
        tt = -rangeT;
        Drift(electron, tt);
        for (int i  = 0; i < nSteps; i++) {
            xsec = xsec_cal(electron, k_vector, wl, tt);
            r = RandFlat->shoot(0.0, 1.0);
            if (r < xsec) {
                sampleScatteredPhoton(electron, k_vector, wl);

                // Determine the properties of the electron after scattering
                const G4double Ei = electron.E;
                electron.E += (ki_fV.x0 - kf_fV.x0); // conservation of energy
                electron.beta = (electron.beta * Ei + ki_fV.x - kf_fV.x) / electron.E; // conservation of momentum
                electron.gamma = electron.E * reciprocal_electron_mass_c2;

                const auto& [posX, posY, posZ] = electron.position + pos0;
                fParticleGun->SetParticlePosition(G4ThreeVector(posX, posY, posZ));
                fParticleGun->GeneratePrimaryVertex(anEvent);
            }
            Drift(electron, dT);
            tt += dT;
        }
    }
}