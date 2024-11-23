#include "Garfield/AvalancheGrid.hh"

#include <TF1.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <numeric>

#include "Garfield/Medium.hh"
#include "Garfield/Random.hh"

namespace Garfield {

void AvalancheGrid::SetGrid(const double xmin, const double xmax,
                            const int xsteps, const double ymin,
                            const double ymax, const int ysteps,
                            const double zmin, const double zmax,
                            const int zsteps) {
  m_gridset = true;

  if (zmin >= zmax || zsteps <= 0 || xmin > xmax || xsteps <= 0 ||
      ymin > ymax || ysteps <= 0) {
    std::cerr << m_className
              << "::SetGrid: Error. Grid is not properly defined.\n";
    return;
  }

  // Creating the z-coordinate grid.
  m_zsteps = zsteps;
  m_zStepSize = (zmax - zmin) / zsteps;
  m_zgrid.clear();
  for (int i = 0; i < zsteps; i++) {
    m_zgrid.push_back(zmin + i * m_zStepSize);
  }
  // Idem to for the y-coordinate grid.
  m_ysteps = ysteps;
  m_yStepSize = (ymax - ymin) / ysteps;
  if (m_yStepSize == 0) m_yStepSize = 1;
  m_ygrid.clear();
  for (int i = 0; i < ysteps; i++) {
    m_ygrid.push_back(ymin + i * m_yStepSize);
  }
  // Idem for the x-coordinate grid.
  m_xsteps = xsteps;
  m_xStepSize = (xmax - xmin) / xsteps;
  if (m_xStepSize == 0) m_xStepSize = 1;
  m_xgrid.clear();
  for (int i = 0; i < xsteps; i++) {
    m_xgrid.push_back(xmin + i * m_xStepSize);
  }

  if (m_debug) {
    std::cout << m_className << "::SetGrid: Grid created:\n"
              << "       x range = (" << xmin << "," << xmax << ").\n"
              << "       y range = (" << ymin << "," << ymax << ").\n"
              << "       z range = (" << zmin << "," << zmax << ").\n";
  }
}

int AvalancheGrid::GetAvalancheSize(double dx, const int nsize,
                                    const double alpha, const double eta) {
  // Algorithm to get the size of the avalanche after it has propagated over a
  // distance dx.

  int newnsize = 0;
  const double k = eta / alpha;
  const double ndx = exp((alpha - eta) * dx);  
  // If the size is higher than 1000 the central limit theorem will be used 
  // to describe the growth of the Townsend avalanche.
  if (nsize < 1000) {
    // Condition to which the random number will be compared. 
    // If the number is smaller than the condition, nothing happens. 
    // Otherwise, the single electron will be attached or retrieve 
    // additional electrons from the gas.
    const double prob = k * (ndx - 1) / (ndx - k);
    // Running over all electrons in the avalanche.
    for (int i = 0; i < nsize; i++) {
      // Draw a random number from the uniform distribution (0,1).
      double s = RndmUniformPos();
      if (s >= prob) {
        newnsize += (int)(1 + log((ndx - k) * (1 - s) / (ndx * (1 - k))) /
                              log(1 - (1 - k) / (ndx - k)));
      }
    }
  } else {
    // Central limit theorem.
    const double sigma = sqrt((1 + k) * nsize * ndx * (ndx - 1) / (1 - k));
    newnsize = RndmGaussian(nsize * ndx, sigma);
  }
  return newnsize;
}

bool AvalancheGrid::SnapToGrid(const double x, const double y,
                               const double z, const double /*v*/,
                               const int n) {
  // Snap electron from AvalancheMicroscopic to the predefined grid.
  if (!m_gridset) {
    std::cerr << m_className << "::SnapToGrid: Error. Grid is not defined.\n";
    return false;
  }
  // Finding the position on the grid.
  // TODO: Snap must be dependent on the direction of drift.
  int indexX = round((x - m_xgrid.front()) / m_xStepSize);
  int indexY = floor((y - m_ygrid.front()) / m_yStepSize);
  int indexZ = round((z - m_zgrid.front()) / m_zStepSize);

  if (m_debug) {
    std::cout << m_className << "::SnapToGrid: ix = " << indexX
              << ", iy = " << indexY << ", iz = " << indexZ << ".\n\n";
  }
  if (indexX < 0 || indexX >= m_xsteps || indexY < 0 || indexY >= m_ysteps ||
      indexZ < 0 || indexZ >= m_zsteps) {
    if (m_debug)
      std::cerr << m_className << "::SnapToGrid: Point is outside the grid.\n";
    return false;
  }

  AvalancheNode newNode;
  newNode.ix = indexX;
  newNode.iy = indexY;
  newNode.iz = indexZ;
  if (!GetParameters(newNode)) {
    // TODO:Memory leak?
    if (m_debug)
      std::cerr << m_className
                << "::SnapToGrid: Could not retrieve parameters from sensor.\n";
    return false;
  }

  // When snapping the electron to the grid the distance traveled can yield
  // additional electrons or attachment.

  double step = z - m_zgrid[indexZ];
  if (newNode.velNormal[0] != 0) {
    step = x - m_xgrid[indexX];
  } else if (newNode.velNormal[1] != 0) {
    step = y - m_ygrid[indexY];
  }

  const int nholder =
      GetAvalancheSize(step, n, newNode.townsend, newNode.attachment);
  if (nholder == 0) {
    if (m_debug)
      std::cerr << m_className << "::SnapToGrid: n from 1 to 0 -> cancel.\n";
    return false;
  }

  newNode.n = nholder < m_MaxSize ? nholder : m_MaxSize;
  m_nTotal += newNode.n;

  bool alreadyExists = false;

  for (AvalancheNode &existingNode : m_activeNodes) {
    if (existingNode.ix == newNode.ix && existingNode.iy == newNode.iy &&
        existingNode.iz == newNode.iz) {
      alreadyExists = true;
      existingNode.n += newNode.n;
    }
  }

  // TODO: What if time changes as you are importing avalanches?
  newNode.time = m_time;
  if (!alreadyExists) m_activeNodes.push_back(newNode);

  if (m_debug) {
    std::cerr << m_className << "::SnapToGrid: n from 1 to " << nholder << ".\n"
              << "    Snapped to (x,y,z) = (" << x << " -> " << m_xgrid[indexX]
              << ", " << y << " -> " << m_ygrid[indexY] << ", " << z << " -> "
              << m_zgrid[indexZ] << ").\n";
  }
  return true;
}

void AvalancheGrid::NextAvalancheGridPoint() {
  // This main function propagates the electrons and applies the avalanche
  // statistics.
  int Nholder = 0;  // Holds the avalanche size before propagating it to the
  // next point in the grid.
  m_run = false;
  for (AvalancheNode &node : m_activeNodes) {  // For every avalanche node

    if (!node.active) {
      continue;
    } else {
      m_run = true;
    }

    if (m_debug)
      std::cerr << m_className << "::NextAvalancheGridPoint:(ix,iy,iz) = ("
                << node.ix << "," << node.iy << "," << node.iz << ").\n";

    // Get avalanche size.
    Nholder = node.n;
     if(node.path.xs.empty()){
          node.path.xs.push_back({m_xgrid[node.ix],
                                    m_ygrid[node.iy],
                                    m_zgrid[node.iz]});
          node.path.ts.push_back(node.time + node.dt);
          node.path.qs.push_back((Nholder) / 2);
      }

    if (Nholder == 0) continue;  // If empty go to next point.
    // If the total avalanche size is smaller than the set saturation
    // limit the GetAvalancheSize function is utilized to obtain the size
    // after its propagation to the next z-coordinate grid point. Else,
    // the size will be kept constant under the propagation.

    if (!m_layerIndix && m_nTotal < m_MaxSize) {
      int holdnsize = GetAvalancheSize(node.stepSize, node.n, node.townsend,
                                       node.attachment);

      if (m_MaxSize - m_nTotal < holdnsize - node.n)
        holdnsize = m_MaxSize - m_nTotal + node.n;

      node.n = holdnsize;

    } else if (m_layerIndix && m_NLayer[node.layer - 1] < m_MaxSize) {
      int holdnsize = GetAvalancheSize(node.stepSize, node.n, node.townsend,
                                       node.attachment);

      if (m_MaxSize - m_NLayer[node.layer - 1] < holdnsize - node.n)
        holdnsize = m_MaxSize - m_NLayer[node.layer - 1] + node.n;

      node.n = holdnsize;

    } else {
      m_Saturated = true;

      if (m_SaturationTime == -1) m_SaturationTime = node.time + node.dt;
    }
    // Produce induced signal on readout electrodes.
    
    node.path.xs.push_back({m_xgrid[node.ix + node.velNormal[0]],
                              m_ygrid[node.iy + node.velNormal[1]],
                              m_zgrid[node.iz + node.velNormal[2]]});
    node.path.ts.push_back(node.time + node.dt);
    node.path.qs.push_back((Nholder + node.n) / 2);
      
    // Update total number of electrons.

    if (m_layerIndix) m_NLayer[node.layer - 1] += node.n - Nholder;

    m_nTotal += node.n - Nholder;

    if (m_diffusion) {
      // TODO: to implement
    }

    if (m_debug) std::cerr << "n = " << Nholder << " -> " << node.n << ".\n";

    // Update position index.

    node.ix += node.velNormal[0];
    node.iy += node.velNormal[1];
    node.iz += node.velNormal[2];

    // After all active grid points have propagated, update the time.
    if (m_debug) std::cerr << "t = " << node.time << " -> ";
    node.time += node.dt;
    if (m_debug) std::cerr << node.time << ".\n";

    DeactivateNode(node);
  }
  if (m_debug) std::cerr << "N = " << m_nTotal << ".\n\n";
}

void AvalancheGrid::DeactivateNode(AvalancheNode &node) {

  if (node.n == 0) node.active = false;

  if (node.velNormal[2] != 0) {
    if ((node.velNormal[2] < 0 && node.iz == 0) ||
        (node.velNormal[2] > 0 && node.iz == m_zsteps - 1))
      node.active = false;
  } else if (node.velNormal[1] != 0) {
    if ((node.velNormal[1] < 0 && node.iy == 0) ||
        (node.velNormal[1] > 0 && node.iy == m_ysteps - 1))
      node.active = false;
  } else {
    if ((node.velNormal[0] < 0 && node.ix == 0) ||
        (node.velNormal[0] > 0 && node.ix == m_xsteps - 1))
      node.active = false;
  }

  double e[3], v;
  int status;
  Medium *m = nullptr;
  m_sensor->ElectricField(m_xgrid[node.ix], m_ygrid[node.iy],
                          m_zgrid[node.iz], e[0], e[1], e[2], v, m,
                          status);

  if (status == -5 || status == -6) {
    node.active = false;  // If not inside a gas gap return false to terminate
  }
        
  if(!node.active){
    // If node has terminated then the signal from the avalanche is calculated
    m_sensor->AddSignalWeightingPotential(-1,node.path.ts,node.path.xs,
                                          node.path.qs);
    if (m_debug)
        std::cerr << m_className << "::DeactivateNode: Node deactivated.\n";
  }
}

void AvalancheGrid::StartGridAvalanche() {
  // Start the AvalancheGrid algorithm.
  if ((!m_importAvalanche && !m_driftAvalanche) || !m_sensor) return;

  std::cerr << m_className
            << "::StartGridAvalanche::Starting grid based simulation with "
            << m_nTotal << " initial electrons.\n";
  if (m_nTotal <= 0) {
    std::cerr << m_className << "::StartGridAvalanche::Cancelled.\n";
    return;
  }

  m_nestart = m_nTotal;

  // Main loop.
  while (m_run == true) {
    if (m_debug)
      std::cout
          << "============ \n" << m_className
          << "::StartGridAvalanche: Looping over nodes.\n ============ \n";
    NextAvalancheGridPoint();
  }

  std::vector<double> tlist = {};
  for (AvalancheNode &node : m_activeNodes) {
    tlist.push_back(node.time);
  }
  double maxTime = *max_element(std::begin(tlist), std::end(tlist));

  if (m_Saturated)
    std::cerr << m_className
              << "::StartGridAvalanche::Avalanche maximum size of " << m_MaxSize
              << " electrons reached at " << m_SaturationTime << " ns.\n";

  std::cerr << m_className
            << "::StartGridAvalanche::Final avalanche size = " << m_nTotal
            << " ended at t = " << maxTime << " ns.\n";

  return;
}

void AvalancheGrid::AvalancheElectron(const double x, const double y,
                                      const double z, const double t,
                                      const int n) {
  m_driftAvalanche = true;

  if (m_time == 0 && m_time != t && m_debug)
    std::cerr << m_className
              << "::CreateAvalanche::Overwriting start time of avalanche for t "
                 "= 0 to " << t << ".\n";
  m_time = t;

  if (SnapToGrid(x, y, z, 0, n) && m_debug)
    std::cerr << m_className
              << "::CreateAvalanche::Electron added at (t,x,y,z) =  (" << t
              << "," << x << "," << y << "," << z << ").\n";
}

void AvalancheGrid::ImportElectronsFromAvalancheMicroscopic(
    AvalancheMicroscopic *avmc) {
  // Get the information of the electrons from the AvalancheMicroscopic class.
  if (!avmc) return;

  if (!m_importAvalanche) m_importAvalanche = true;

  // Get initial positions of electrons
  for (const auto& electron : avmc->GetElectrons()) {
    // If the electron is not stopped due to
    // the upper bound of the time range: then skip this electron.
    if (electron.status != -17) continue;
    const auto& p1 = electron.path.at(0);
    const auto& p2 = electron.path.back();
    const double vel = (p2.z - p1.z) / (p2.t - p1.t);
    m_time = p2.t;

    if (SnapToGrid(p2.x, p2.y, p2.z, vel) && m_debug) {
      std::cerr << m_className
                << "::ImportElectronsFromAvalancheMicroscopic: Electron added at "
                   "(x,y,z) =  (" << p2.x << "," << p2.y << "," << p2.z << ").\n";
    }
  }
}

bool AvalancheGrid::GetParameters(AvalancheNode &node) {
  if (!m_sensor) return false;

  double x = m_xgrid[node.ix];
  double y = m_ygrid[node.iy];
  double z = m_zgrid[node.iz];

  if (m_debug)
    std::cerr << m_className
              << "::GetParametersFromSensor::Getting parameters from "
                 "(x,y,z) =  (" << x << "," << y << "," << z << ").\n";

  double e[3], v;
  int status;
  Medium *m = nullptr;
  m_sensor->ElectricField(x, y, z, e[0], e[1], e[2], v, m, status);

  if (m_debug)
    std::cerr << m_className << "::GetParametersFromSensor::status = " << status
              << ".\n";

  if (status == -5 || status == -6)
    return false;  // If not inside a gas gap return false to terminate

  if (m_Townsend >=
      0) {  // If Townsend coef. is not set by user, take it from the sensor.
    node.townsend = m_Townsend;
  } else {
    m->ElectronTownsend(e[0], e[1], e[2], 0., 0., 0., node.townsend);
  }

  if (m_Attachment >=
      0) {  // If attachment coef. is not set by user, take it from the sensor.
    node.attachment = m_Attachment;
  } else {
    m->ElectronAttachment(e[0], e[1], e[2], 0., 0., 0., node.attachment);
  }

  if (m_Velocity >
      0) {  // If velocity is not set by user, take it from the sensor.
    node.velocity = m_Velocity;
    node.velNormal = m_velNormal;
  } else {
    double vx, vy, vz;
    m->ElectronVelocity(e[0], e[1], e[2], 0., 0., 0., vx, vy, vz);

    double vel = sqrt(vx * vx + vy * vy + vz * vz);
    if (vel == 0.) return false;
    if (vel != std::abs(vx) && vel != std::abs(vy) && vel != std::abs(vz))
      return false;
    int nx = (int)round(vx / vel);
    int ny = (int)round(vy / vel);
    int nz = (int)round(vz / vel);

    node.velNormal = {nx, ny, nz};
    node.velocity = -std::abs(vel);
  }

  if (node.velNormal[0] != 0) {
    node.stepSize = m_xStepSize;
  } else if (node.velNormal[1] != 0) {
    node.stepSize = m_yStepSize;
  } else {
    node.stepSize = m_zStepSize;
  }

  if (m_debug) {
    std::cerr << m_className << "::GetParametersFromSensor:\n"
              << "    stepSize = " << node.stepSize << " [cm].\n"
              << "    velNormal = (" << node.velNormal[0] << ", "
              << node.velNormal[1] << ", " << node.velNormal[2] << ") [1].\n";
  }
  node.dt = std::abs(node.stepSize / node.velocity);

  // print
  if (m_debug || !m_printPar) {
    std::cerr << m_className << "::GetParametersFromSensor:\n"
              << "    Electric field = (" << 1.e-3 * e[0] << ", "
              << 1.e-3 * e[1] << ", " << 1.e-3 * e[2] << ") [kV/cm].\n"
              << "  Townsend = " << node.townsend
              << " [1/cm], Attachment = " << node.attachment
              << " [1/cm], Velocity = " << node.velocity << " [cm/ns].\n";
  }
  if (m_debug)
    std::cerr << m_className << "::StartGridAvalanche::Time steps per loop "
              << node.dt << " ns.\n";
  m_printPar = true;
  return true;
}

void AvalancheGrid::Reset() {

  std::cerr << m_className << "::Reset::Resetting AvalancheGrid.\n";
  m_time = 0;
  m_nTotal = 0;
  m_run = true;

  m_Saturated = false;
  m_SaturationTime = -1;

  m_driftAvalanche = false;

  m_activeNodes.clear();
  m_layerIndix = false;
  m_NLayer.clear();
}

void AvalancheGrid::AsignLayerIndex(ComponentParallelPlate *RPC) {
  int im = 0;
  double epsM = 0;
  std::vector<double> nLayer(RPC->NumberOfLayers(), 0);
  for (AvalancheNode &node : m_activeNodes) {
    double y = m_ygrid[node.iy];
    RPC->getLayer(y, im, epsM);
    node.layer = im;
    nLayer[im - 1] += node.n;
    // std::cerr << m_className << "::AsignLayerIndex::im = "<<im<<".\n";
  }

  m_NLayer = nLayer;
  m_layerIndix = true;
}

}  // namespace Garfield
