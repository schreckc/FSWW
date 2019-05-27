#ifndef WATERSURFACE_HPP
#define WATERSURFACE_HPP

#include <list>

#include "definitions.hpp"
#include "Wave.hpp"
#include "Obstacle.hpp"
#include "Grid.hpp"
#include "CudaWaterSurface.hpp"
#include "ProjectedGrid.hpp"
#include "MovingEquivalentSource.hpp"
#include <thread>

class WaterSurface {
public:
  WaterSurface();
  ~WaterSurface();

  void clear();
  void reset();
  
  void setLists();
  void setAmpli();
  void setAmpli(FLOAT t);
  void setObstacleGrid();

  FLOAT height(uint i, uint j) const;
  FLOAT height_obs(uint i, uint j) const;

  void addEqSource(FLOAT x, FLOAT y, FLOAT wl = 0, COMPLEX ampli = COMPLEX(1, 0));
  void addLinearWave(FLOAT x, FLOAT y, FLOAT wl = 0, COMPLEX ampli = COMPLEX(1, 0));

  void applyForce(FLOAT x, FLOAT y, FLOAT radius);
  void addDrop(FLOAT x, FLOAT y, FLOAT size);
  int getTDrop(FLOAT size, int t);
  void addDropW(VEC2 pos, FLOAT size, int t);
  void update();
  friend void updateSourcesAmplis(WaterSurface *ws);
  void updateHeight();

  void getObstacleBoundary(std::list<VEC2> &boundaries) const;

  void exportAmplitude(std::string file) const;
  void exportAmplitudeIm(std::string file) const;
  void exportAmplitudeRe(std::string file) const;
  void exportPhase(std::string file) const;

  void exportAmplitudeIn(std::string file) const;
  void exportAmplitudeInIm(std::string file) const;
  void exportAmplitudeInRe(std::string file) const;
  void exportPhaseIn(std::string file) const;

  void exportAmplitudeScattered(std::string file) const;
  void exportAmplitudeScatteredIm(std::string file) const;
  void exportAmplitudeScatteredRe(std::string file) const;
  void exportPhaseScattered(std::string file) const;

  void exportDirectivityAmplitude(std::string file, FLOAT dist) const;
  void exportDirectivityAnalytic(std::string file, FLOAT r_obs, FLOAT r) const;

  void exportMitsuba(std::string file) const;
  void exportSurfaceFreq(std::string file) const;
  void importSurfaceFreq(std::string file);
  void exportSurfaceTime(std::string file) const;
  void importSurfaceTime(std::string file);

  void exportSurfaceObj(std::string file) const;
  void exportObstacleGrid(std::string file) const;
  std::ofstream &exportObstacleMitsuba(std::ofstream &file) const;
  
  void importConfig(std::string file);

  void setImport(std::string file);
  void setExport(std::string file);
  void setExportMitsuba(std::string file);
  void setImportConf(std::string file);
  void setExportStep(uint es);
  void setStopTime(uint end);
  void setData(std::string file);
  void drawHeighField(std::string file);
#ifdef USE_CUDA
  void setCuda();
  void updateCuda();
  void setCudaM();
  void updateCudaM();
#endif

#ifdef PROJECTED_GRID
  void updatePosGridCuda();  
  void setProjGrid(uint i, uint j, FLOAT x, FLOAT y);
  VEC3 getPosProjGrid(uint i, uint j) const;
  VEC3 getPosProjGrid(uint i) const;

  void setTargetLookAt(VEC2 target);
#endif
  VEC3 getPosGrid(uint i, uint j) const;
  VEC3 getPosGrid(uint i) const;

  FLOAT minWL() const;
  FLOAT maxWL() const;

  void setObsGridCuda();

  void addMovingSource(VEC2 pos0, VEC2 speed, FLOAT start, FLOAT stop);
private:
  FLOAT step_wl;
  FLOAT min_wl, max_wl;
  uint nb_wl;
  
  Grid u, u_re, u_im;
  Grid b;
  Grid dispx, dispy;
  std::vector<Grid> ampli_in_re;
  std::vector<Grid> ampli_in_im;
  std::vector<Grid> ampli_scattered_re;
  std::vector<Grid> ampli_scattered_im;

  int time;

  std::vector<std::list<Wave*> > waves;
  std::vector<std::list<MovingEquivalentSource*> > moving_waves;
  std::vector<FLOAT> wave_lenghts;

  std::list<Obstacle*> obstacles;

  bool import_;
  bool export_;
  bool export_mit;
  bool load_conf;
  
  std::string import_file;
  std::string export_file;
  std::string export_mit_file;
  std::string conf_file;
  std::string data_file;

  uint stop_time;
  uint export_step;

#ifdef USE_CUDA
  mutable CudaWaterSurface cuda_surface;
#endif
  
  std::vector<int> index_inter_src;
  std::vector<std::vector<EquivalentSource*> > inter_src;

  ProjectedGrid proj_grid;
  std::thread t_solve;

  VEC2 target_lookat;
};



#endif
