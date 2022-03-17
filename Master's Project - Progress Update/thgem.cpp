//*****************************************************
//	Code written by:
//			 Danush Shekar,
//			 NISER, Bhubaneswar, India,
//			 (11th March 2022)
//*****************************************************
#include <iostream>
#include <TApplication.h>
#include <fstream>
#include "Garfield/SolidBox.hh" 
#include "Garfield/SolidHole.hh" 
#include "Garfield/SolidTube.hh" 
#include "Garfield/SolidRidge.hh"
#include "Garfield/GeometrySimple.hh"
#include "Garfield/MediumMagboltz.hh"
#include "Garfield/MediumConductor.hh"
#include "Garfield/MediumPlastic.hh"
#include "Garfield/ComponentNeBem3d.hh"
#include "Garfield/ViewGeometry.hh"
#include "Garfield/ViewField.hh"
#include "Garfield/Random.hh"
#include "Garfield/Sensor.hh"
#include "Garfield/AvalancheMC.hh"
#include "Garfield/AvalancheMicroscopic.hh"
#include "Garfield/TrackHeed.hh"
#include "Garfield/ViewDrift.hh"
#include "Garfield/ViewSignal.hh"

using namespace std;
using namespace Garfield;

int main(int argc, char * argv[]) {

  TApplication app("app", &argc, argv);

// Define materials
// Gas mixture
  MediumMagboltz *gas = new MediumMagboltz();
  const double pressure = 760.;
  const double temperature = 293.15;
  gas->SetTemperature(temperature);
  gas->SetPressure(pressure);
  gas->SetComposition("Ar", 90., "CO2", 10.);
  gas->LoadGasFile("../ar_90_co2_10_80k.gas");


// Read the ion mobility table from file.
  const string garfpath = getenv("GARFIELD_HOME");
  gas->LoadIonMobility(garfpath + "/Data/IonMobility_Ar+_Ar.txt");

// Other mediums
  MediumConductor Cu;
  MediumPlastic kp;
  kp.SetDielectricConstant(3.9);

// Geometry.
  GeometrySimple geo;
  geo.SetMedium(gas);

  double i_dia = 0.0200;
  double o_dia = 0.0300;
  double pitch = 0.0450;
  double ylen = sqrt(3)*pitch/2, xlen = pitch;
  //double vpitch = sqrt(3)*pitch/2, hpitch = pitch/2;
  double cu_height = 0.0005;
  double kp_height = 0.0250;
  double drift_gap = 0.5;
  double ind_gap = 0.25;

  double anode_z = 0.0 + cu_height/2;
  double drift_z = anode_z + ind_gap + 2*cu_height + kp_height + drift_gap + cu_height;
  double cu_upper_z = anode_z + ind_gap + cu_height + kp_height + cu_height;
  double kapton_z = anode_z + ind_gap + cu_height + kp_height/2 + cu_height/2;
  double cu_lower_z = anode_z + ind_gap + cu_height;

  double position_x[3] = {-pitch, 0.0, pitch};
  double position_y[3] = {-pitch, 0.0, pitch};

  double driftV = -2500;
  double upperV = -1500;
  double lowerV = -500;
  double anodeV = 0;
  
  double offset_x = 0.0, offset_y = 0.0, offset_z = 0.0;//cu_height/2 + kp_height/2;
  double offset_x1 = -xlen/4, offset_y1 = -ylen/2, offset_z1 = 0.0;
  double offset_x2 = +xlen/4, offset_y2 = ylen/2, offset_z2 = 0.0;
  //double offset_x3 = +hpitch/2, offset_y3 = +vpitch/2, offset_z3 = 0.0;
  //double offset_x4 = -hpitch/2, offset_y4 = +vpitch/2, offset_z4 = 0.0;
  
  std::cout<<"\nPitch = "<<pitch<<endl;
  
  SolidBox drift(offset_x + offset_x1, offset_y + offset_y1, drift_z, xlen/2, ylen/2, cu_height/2);
  drift.SetBoundaryPotential(driftV);
  geo.AddSolid(&drift, &Cu);
  
  SolidBox drift2(offset_x + offset_x2, offset_y + offset_y2, drift_z, xlen/2, ylen/2, cu_height/2);
  drift2.SetBoundaryPotential(driftV);
  geo.AddSolid(&drift2, &Cu);
  
  
  
  SolidHole cu_upper(offset_x + offset_x1, offset_y + offset_y1, cu_upper_z, o_dia/2, o_dia/2, xlen/2, ylen/2, cu_height/2);
  cu_upper.SetBoundaryPotential(upperV);
  geo.AddSolid(&cu_upper, &Cu);

  SolidHole kapton(offset_x + offset_x1, offset_y + offset_y1, kapton_z, i_dia/2, i_dia/2, xlen/2, ylen/2, kp_height/2);
  kapton.SetBoundaryDielectric();
  geo.AddSolid(&kapton, &kp);
  
  SolidHole cu_lower(offset_x + offset_x1, offset_y + offset_y1, cu_lower_z, o_dia/2, o_dia/2, xlen/2, ylen/2, cu_height/2);
  cu_lower.SetBoundaryPotential(lowerV);
  geo.AddSolid(&cu_lower, &Cu);

 
  
  SolidHole cu_upper2(offset_x + offset_x2, offset_y + offset_y2, cu_upper_z, o_dia/2, o_dia/2, xlen/2, ylen/2, cu_height/2);
  cu_upper2.SetBoundaryPotential(upperV);
  geo.AddSolid(&cu_upper2, &Cu);

  SolidHole kapton2(offset_x + offset_x2, offset_y + offset_y2, kapton_z, i_dia/2, i_dia/2, xlen/2, ylen/2, kp_height/2);
  kapton2.SetBoundaryDielectric();
  geo.AddSolid(&kapton2, &kp);
  
  SolidHole cu_lower2(offset_x + offset_x2, offset_y + offset_y2, cu_lower_z, o_dia/2, o_dia/2, xlen/2, ylen/2, cu_height/2);
  cu_lower2.SetBoundaryPotential(lowerV);
  geo.AddSolid(&cu_lower2, &Cu);
  
  

  SolidBox Anode(offset_x + offset_x1, offset_y + offset_y1, offset_z + anode_z, xlen/2, ylen/2, cu_height/2);
  Anode.SetBoundaryPotential(anodeV);
  Anode.SetLabel("anode");
  geo.AddSolid(&Anode, &Cu);
  
  SolidBox Anode2(offset_x + offset_x2, offset_y + offset_y2, offset_z + anode_z, xlen/2, ylen/2, cu_height/2);
  Anode2.SetBoundaryPotential(anodeV);
  Anode2.SetLabel("anode2");
  geo.AddSolid(&Anode2, &Cu);
  
  
  //=============================================================================================================================================
  
  
  // Plot device geometry in 3D
  /*ViewGeometry geomView;
  geomView.SetGeometry(&geo);
  geomView.Plot();
  app.Run();	// Keep this on to interact with the 3D figure
  */

  // Plot device geometry in 2D
  /*
  ViewGeometry geomView2d;
  geomView2d.SetGeometry(&geo);
  geomView2d.SetArea(-2*pitch, -2*pitch, -10*kp_height, 2*pitch, 2*pitch, 5*kp_height);
  geomView2d.SetPlane(0, 1, 0, 0, 0, 0.0);
  geomView2d.Plot2d();
  */
  
  double tgtElSize = 1.e-2;	// target element size
  int minEl = 3, maxEl = 8;	// minimum and maximum number of elements
  int xcopy = 40, ycopy = 40, zcopy = 0;	// no. of copies in the 3 directions
  ComponentNeBem3d nebem;
  nebem.SetGeometry(&geo);
  nebem.SetNumberOfThreads(3);	// Set no. of threads for the calculation
  nebem.SetTargetElementSize(tgtElSize);
  nebem.SetMinMaxNumberOfElements(minEl, maxEl);
  nebem.SetPeriodicityX(xlen);
  nebem.SetPeriodicityY(2*ylen);
  nebem.SetPeriodicCopies(xcopy, ycopy, zcopy);
  nebem.UseLUInversion();
  //nebem.EnableDebugging();
  nebem.Initialise();

  int plot_field = 1;
  int calc_field = 0;
  
  if(plot_field==1)
  {
    std::cout<<"\n Plotting of equipotential lines calculations have begun."<<endl;
    ViewField fieldView;
    fieldView.SetComponent(&nebem);
    fieldView.SetNumberOfContours(100);
    fieldView.SetPlane(0, -1, 0, offset_x2, offset_y2, 0);     // Set the normal vector of the viewing plane (xz plane).
    fieldView.SetArea(-1*xlen, cu_lower_z - 0.1, 1*xlen, cu_upper_z + 0.1); // Set the plot limits in the current viewing plane.
    TCanvas* cf = new TCanvas("cf", "Potential (V)", 600, 500);
    //cf->SetLeftMargin(0.16);
    fieldView.SetCanvas(cf);
    fieldView.PlotContour("v");
    char name[1024];
    sprintf(name, "Potential_full(%dV, %dV, %dV, %dV).pdf", (int)anodeV, (int)lowerV, (int)upperV, (int)driftV);
    cf->SaveAs(name);
  }
  
  if(calc_field==1){
    std::cout<<"\n Field calculations have begun."<<endl;
    {// field along line 1
    std::ofstream fldfile;
    fldfile.open("40x40_unit_cell_fieldline_along_x(y=0.5drift_gap).csv");	
    
    int nx = 1000;
    double delx = ((2*xcopy+1)*xlen) / (double)(nx - 1);
    double xp, yp = offset_y2, zp = offset_z + cu_upper_z + drift_gap/2;
    Medium* medium = nullptr; 
    double ex = 0., ey = 0., ez = 0., e = 0.0, v = 0.;
    int status = 0;
    for(int ix = 0; ix < nx; ++ix)
    {
      xp = offset_x - xcopy*xlen + ix*delx;
      nebem.ElectricField(xp, yp, zp, ex, ey, ez, v, medium, status);
      e = (ex*ex + ey*ey + ez*ez);
      e = pow(e, 0.5);
      fldfile << xp << "," << yp << "," << zp << "," << ex << "," << ey << "," << ez << "," << e << "," << v << "," << medium << "," << status << std::endl;
		}
    fldfile.close();
    }
    
    {// field along line 2
    std::ofstream fldfile;
    fldfile.open("40x40_unit_cell_fieldline_along_x(y=0.05_from_top_foil).csv");	
    
    int nx = 1000;
    double delx = ((2*xcopy+1)*xlen) / (double)(nx - 1);
    double xp, yp = offset_y2, zp = offset_z + cu_upper_z + drift_gap/10;
    Medium* medium = nullptr; 
    double ex = 0., ey = 0., ez = 0., e = 0.0, v = 0.;
    int status = 0;
    for(int ix = 0; ix < nx; ++ix)
    {
      xp = offset_x - xcopy*xlen + ix*delx;
      nebem.ElectricField(xp, yp, zp, ex, ey, ez, v, medium, status);
      e = (ex*ex + ey*ey + ez*ez);
      e = pow(e, 0.5);
      fldfile << xp << "," << yp << "," << zp << "," << ex << "," << ey << "," << ez << "," << e << "," << v << "," << medium << "," << status << std::endl;
		}
    fldfile.close();
    }
    
    
    {// field along line 3
    std::ofstream fldfile;
    fldfile.open("40x40_unit_cell_fieldline_along_x(y=0.5ind_gap).csv");	
    
    int nx = 1000;
    double delx = ((2*xcopy+1)*xlen) / (double)(nx - 1);
    double xp, yp = offset_y2, zp = cu_lower_z - ind_gap/2;
    Medium* medium = nullptr; 
    double ex = 0., ey = 0., ez = 0., e = 0.0, v = 0.;
    int status = 0;
    for(int ix = 0; ix < nx; ++ix)
    {
      xp = offset_x - xcopy*xlen + ix*delx;
      nebem.ElectricField(xp, yp, zp, ex, ey, ez, v, medium, status);
      e = (ex*ex + ey*ey + ez*ez);
      e = pow(e, 0.5);
      fldfile << xp << "," << yp << "," << zp << "," << ex << "," << ey << "," << ez << "," << e << "," << v << "," << medium << "," << status << std::endl;
		}
    fldfile.close();
    }
    
    {// field along line 4
    std::ofstream fldfile;
    fldfile.open("40x40_unit_cell_fieldline_along_x(y=0.05_from_bottom_foil).csv");	
    
    int nx = 1000;
    double delx = ((2*xcopy+1)*xlen) / (double)(nx - 1);
    double xp, yp = offset_y2, zp = cu_lower_z - ind_gap/4;
    Medium* medium = nullptr; 
    double ex = 0., ey = 0., ez = 0., e = 0.0, v = 0.;
    int status = 0;
    for(int ix = 0; ix < nx; ++ix)
    {
      xp = offset_x - xcopy*xlen + ix*delx;
      nebem.ElectricField(xp, yp, zp, ex, ey, ez, v, medium, status);
      e = (ex*ex + ey*ey + ez*ez);
      e = pow(e, 0.5);
      fldfile << xp << "," << yp << "," << zp << "," << ex << "," << ey << "," << ez << "," << e << "," << v << "," << medium << "," << status << std::endl;
		}
    fldfile.close();
    }

  cout<<"\n\n :: Complete :: \n\n";
  app.Run();
  return 0;
  } 
}
