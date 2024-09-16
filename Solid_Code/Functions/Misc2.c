
inline void ReportData()
{
  if (times.use_real_time_to_save)
    {
      if(times.real_time >= times.time_print) 
      {
        ReportStatus(); 
      }
      if(times.real_time >= times.time_save_csv) 
      {
        AppendToCSVDataFile();
        ReportBoundary();
        times.time_save_csv += times.dtime_save_csv;
      }
      if(times.real_time >= times.time_save_vtk) 
      {
        WriteVTKDataFile((int) (times.time_save_vtk/times.dtime_save_vtk));
        times.time_save_vtk += times.dtime_save_vtk;
      }
    }
  else
  {
    if(times.time_step % times.ntimeSteps_print == 0) 
    {
      ReportStatus(); 
    }
    if(times.time_step % times.ntimeSteps_save_csv == 0)
    { 
      AppendToCSVDataFile();
      ReportBoundary();
    }
    if(times.time_step % times.ntimeSteps_save_vtk == 0) 
    {
      WriteVTKDataFile((int) times.time_step/times.ntimeSteps_save_vtk);
    }
  }
}