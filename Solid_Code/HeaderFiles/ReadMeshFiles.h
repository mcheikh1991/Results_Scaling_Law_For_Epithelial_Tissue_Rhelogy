
/* A header files that reads:
    x, y, z coodinates of each point
    line information
    triangle information
    square information
    LE.K_line information
  */

//

// Reading the points
if (vtkInputFile[0] != '\0')  ReadPointsFromVTKFile(vtkInputFile);
else                          ReadPointsFromMeshFolder(meshInputFolder);

if (print_mesh_points)
  for(int i=0; i<points.N; i++)
    printf("%d: (%lf, %lf, %lf)\n",i, points.x[i], points.y[i],points.z[i]);


// Creating force for each point
forces.N = points.N;
forces.x = (double *)calloc(forces.N,sizeof(double));
forces.y = (double *)calloc(forces.N,sizeof(double));
forces.z = (double *)calloc(forces.N,sizeof(double));

//===========================================================================
//===========================================================================

if (vtkInputFile[0] != '\0') ReadCellsFromVTKFile(vtkInputFile);
else                         ReadCellsFromMeshFolder(meshInputFolder);


if (print_mesh_lines) {
  for(int i=0; i<lines.N; i++)      printf("%d: (%d, %d)\n",i, lines.I0[i], lines.I1[i]);
  for(int i=0; i<triangles.N; i++)  printf("%d: (%d, %d, %d)\n",i, triangles.I0[i], triangles.I1[i], triangles.I2[i]);
  for(int i=0; i<squares.N; i++)    printf("%d: (%d, %d, %d, %d)\n",i, squares.I0[i], squares.I1[i], squares.I2[i], squares.I3[i]);
}

//===========================================================================
// Other Mesh Data
//===========================================================================
ReadLinearElasticMeshData();
ReadSurfaceTensionMeshData();
ReadVolumeConservationMeshData();
ReadL0DynamicsMeshData();