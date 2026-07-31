//detector geometry needed in analysis

#ifndef GEOMETRY_H
#define GEOMETRY_H

namespace geometry{
	
	//TPC parameters
	constexpr double tpc_center_x = 0;
	constexpr double tpc_center_y = 0;
	constexpr double tpc_center_z = 0;
	constexpr double tpcInstrumentedRadius = 260; //cm
	constexpr double tpcInstrumentedLength = 500; //cm
	constexpr double B_field = 0.5; //T
	constexpr double adc_sampling_frequency = 20; //MHz
	constexpr double spatial_resolution = 2; //mm, could be 1mm, yz resolution
	constexpr double pixel_spacing_cm = 0.4; //cm, from NDLAr
	constexpr double X0 = 11.93; //m, radiation length of gaseous argon, from Fredrico's Kalman Filter Paper
	constexpr double drift_velocity = 3.011; //of electrons in gas, cm/microsecond, from Naseem's slides
	
	//ECal parameters
	constexpr double ecalLayers = 42;
	
	//MuID parameters
	constexpr double muidLayers = 3;

}

#endif
