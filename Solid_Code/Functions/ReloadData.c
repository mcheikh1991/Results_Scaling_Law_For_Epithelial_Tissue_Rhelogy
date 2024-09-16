
void InitailizeOrReloadData()
{
  if(times.reload) // Reload
  {

    if (reloadVTKFile[0] != '\0')  ReloadPointsFromVTKFile(reloadVTKFile);
    else                           ReloadPointsFromMeshFolder(reloadFolder);

    times.real_time = times.reload_time;
    times.time_save_vtk = times.real_time + times.dtime_save_vtk; 
    times.time_save_csv = times.real_time + times.dtime_save_csv; 
    times.time_print =    times.real_time + times.dtime_print; 
    times.time_step = (int) times.real_time/times.dt;

    if (LD.Type!=L0_DYNAMICS_NULL)    
    {
      if (reloadVTKFile[0] != '\0') ReloadL0FromVTKFile(reloadVTKFile);
      else                          ReloadL0FromMeshFolder(reloadFolder);
    }

  }
  else // Initialize
  {
    CreateCSVDataFile();
    times.time_save_vtk = 0; 
    times.time_save_csv = 0; 
    times.time_print = 0; 
    times.time_step = 0;
  }
}