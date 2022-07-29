/*
###############################################################################
# If you use PhysiCell in your project, please cite PhysiCell and the version #
# number, such as below:                                                      #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1].    #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# See VERSION.txt or call get_PhysiCell_version() to get the current version  #
#     x.y.z. Call display_citations() to get detailed information on all cite-#
#     able software used in your PhysiCell application.                       #
#                                                                             #
# Because PhysiCell extensively uses BioFVM, we suggest you also cite BioFVM  #
#     as below:                                                               #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1],    #
# with BioFVM [2] to solve the transport equations.                           #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# [2] A Ghaffarizadeh, SH Friedman, and P Macklin, BioFVM: an efficient para- #
#     llelized diffusive transport solver for 3-D biological simulations,     #
#     Bioinformatics 32(8): 1256-8, 2016. DOI: 10.1093/bioinformatics/btv730  #
#                                                                             #
###############################################################################
#                                                                             #
# BSD 3-Clause License (see https://opensource.org/licenses/BSD-3-Clause)     #
#                                                                             #
# Copyright (c) 2015-2021, Paul Macklin and the PhysiCell Project             #
# All rights reserved.                                                        #
#                                                                             #
# Redistribution and use in source and binary forms, with or without          #
# modification, are permitted provided that the following conditions are met: #
#                                                                             #
# 1. Redistributions of source code must retain the above copyright notice,   #
# this list of conditions and the following disclaimer.                       #
#                                                                             #
# 2. Redistributions in binary form must reproduce the above copyright        #
# notice, this list of conditions and the following disclaimer in the         #
# documentation and/or other materials provided with the distribution.        #
#                                                                             #
# 3. Neither the name of the copyright holder nor the names of its            #
# contributors may be used to endorse or promote products derived from this   #
# software without specific prior written permission.                         #
#                                                                             #
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" #
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE   #
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE  #
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE   #
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR         #
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF        #
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS    #
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN     #
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)     #
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE  #
# POSSIBILITY OF SUCH DAMAGE.                                                 #
#                                                                             #
###############################################################################
*/

#include "./custom.h"
static double four_thirds_pi =  4.188790204786391;
// ellipse radius
double ellipse_radius= 1.0;
// ellipse matrix
std::vector<std::vector<double>> matrix{{1,0,0},{0,1,0},{0,0,1}};
//
std::vector<double> a_axis(3,0.0);
std::vector<double> b_axis(3,0.0);
void create_cell_types( void )
{
	


	// set the random seed 
	SeedRandom( parameters.ints("random_seed") );  
	
	/* 
	   Put any modifications to default cell definition here if you 
	   want to have "inherited" by other cell types. 
	   
	   This is a good place to set default functions. 
	*/ 
	
	initialize_default_cell_definition(); 
	cell_defaults.phenotype.secretion.sync_to_microenvironment( &microenvironment ); 
	
	cell_defaults.functions.volume_update_function = standard_volume_update_function;
	cell_defaults.functions.update_velocity = standard_update_cell_velocity;

	cell_defaults.functions.update_migration_bias = NULL; 
	cell_defaults.functions.update_phenotype = NULL; // update_cell_and_death_parameters_O2_based; 
	cell_defaults.functions.custom_cell_rule = NULL; 
	cell_defaults.functions.contact_function = NULL; 
	
	cell_defaults.functions.add_cell_basement_membrane_interactions = NULL; 
	cell_defaults.functions.calculate_distance_to_membrane = NULL; 
	
	/*
	   This parses the cell definitions in the XML config file. 
	*/
	
	initialize_cell_definitions_from_pugixml(); 

	/*
	   This builds the map of cell definitions and summarizes the setup. 
	*/
		
	build_cell_definitions_maps(); 

	/*
	   This intializes cell signal and response dictionaries 
	*/

	setup_signal_behavior_dictionaries(); 	

	/* 
	   Put any modifications to individual cell definitions here. 
	   
	   This is a good place to set custom functions. 
	*/ 
	
	cell_defaults.functions.update_phenotype = phenotype_function; 
	cell_defaults.functions.custom_cell_rule = custom_function; 
	cell_defaults.functions.contact_function = contact_function; 
	
	/*
	   This builds the map of cell definitions and summarizes the setup. 
	*/
	// void create_cargo_cell_type( void ) 
	// {
	// 	cargo_cell = find_cell_definition( "cargo cell" ); 

	// 	int apoptosis_index = cell_defaults.phenotype.death.find_death_model_index( PhysiCell_constants::apoptosis_death_model ); 
	
	// 	int oxygen_ID = microenvironment.find_density_index( "oxygen" ); // 0 
	// 	int attract_ID = microenvironment.find_density_index( "chemoattractant" ); // 1
	// 	int therapy_ID = microenvironment.find_density_index( "therapeutic" ); // 2
	
	// 	// reduce o2 uptake 
	// pCell->phenotype.geometry.axis_a
	
	// double axis_a=pCell-> parameters.doubles("cx")
	// double axis_b=pCell-> parameters.doubles("cy")
	// cargo_cell->phenotype.secretion.uptake_rates[oxygen_ID] *= 
	// 	parameters.doubles("cargo_o2_relative_uptake");   
	
	// cargo_cell->phenotype.mechanics.cell_cell_adhesion_strength *= 
	// 	parameters.doubles("cargo_relative_adhesion"); // 0.0;
	// cargo_cell->phenotype.mechanics.cell_cell_repulsion_strength *= 
	// 	parameters.doubles("cargo_relative_repulsion"); // 5.0;
		
	// // figure out mechanics parameters 
	
	// cargo_cell->phenotype.mechanics.relative_maximum_attachment_distance 
	// 	= parameters.doubles( "max_attachment_distance" ) / cargo_cell->phenotype.geometry.radius ; 

	// cargo_cell->phenotype.mechanics.relative_detachment_distance 
	// 	= parameters.doubles( "max_elastic_displacement" ) / cargo_cell->phenotype.geometry.radius ; 
		
	// cargo_cell->phenotype.mechanics.attachment_elastic_constant 
	// 	= parameters.doubles("elastic_coefficient"); 
	
	// // set functions 
	
	// cargo_cell->functions.update_phenotype = cargo_cell_phenotype_rule; 
	// cargo_cell->functions.custom_cell_rule = cargo_cell_rule; 
	// cargo_cell->functions.contact_function = biorobots_contact_function; 
	// cargo_cell->functions.update_migration_bias = NULL;	
	
	// // set custom data values 
	
	// return;
	// }	
	display_cell_definitions( std::cout ); 
	
	return; 
}

void setup_microenvironment( void )
{
	// set domain parameters 
	
	// put any custom code to set non-homogeneous initial conditions or 
	// extra Dirichlet nodes here. 
	
	// initialize BioFVM 
	
	initialize_microenvironment(); 	
	
	return; 
}
double custom_volume_update(double a, double b, double c )
{
	double ellipsoid_volume= four_thirds_pi*a*b*c;
	return ellipsoid_volume;
}
//-----------------------------------------
void setup_tissue( void )
{
	double Xmin = microenvironment.mesh.bounding_box[0]; 
	double Ymin = microenvironment.mesh.bounding_box[1]; 
	double Zmin = microenvironment.mesh.bounding_box[2]; 

	double Xmax = microenvironment.mesh.bounding_box[3]; 
	double Ymax = microenvironment.mesh.bounding_box[4]; 
	double Zmax = microenvironment.mesh.bounding_box[5]; 
	
	if( default_microenvironment_options.simulate_2D == true )
	{
		Zmin = 0.0; 
		Zmax = 0.0; 
	}
	
	double Xrange = Xmax - Xmin; 
	double Yrange = Ymax - Ymin; 
	double Zrange = Zmax - Zmin; 

    double radius = 0.4 * Xrange;
	Cell* pC;
	
    std::vector<double> position = {0,0,0}; 
    static double pi = 3.141592653589793238462643383279502884;
    // double theta = 0.0;
    double theta = 4.2;  // about 2/3 of 2pi
    // double theta_del = 2.0*pi / parameters.ints("number_of_cells");
    double theta_del = -2.0*pi / 3.0 / parameters.ints("number_of_cells");
    position[2] = 0.0;
	for( int k=0; k <parameters.ints("number_of_cells") ; k++ )
	{
		// Cell_Definition* pCD = cell_definitions_by_index[k]; 
		Cell_Definition* pCD = cell_definitions_by_name["ellipsey"]; 

		double semimajor = parameters.doubles("major_axis_2a")/2;
		double ecc = parameters.doubles("eccentricity");
		double vol = parameters.doubles("starting_cell_volume");
		// double b_axis_calc = pow( pow(semimajor,2)*(1-pow(ecc,2)), 0.5);
		// double c_axis_calc = (3*vol)/( 4*pi*pow(((1-ecc)*pow(semimajor,2) ),0.5));
        // Didi update 7-29-22
        double b_axis_calc = semimajor*pow( (1.0-pow(ecc,2)), 0.5);
        double c_axis_calc = (3.0*vol)/( 4.0*pi*pow(semimajor,2)*pow((1-ecc),0.5));

		std::cout << "ecc " << ecc << " ... " << std::endl;
		std::cout << "bax " << b_axis_calc << " ... " << std::endl; 

        // create 4 cells around the origin (in Z=0 plane)
        pC = create_cell( *pCD ); 
        position[0] = radius * std::cos(theta);
        position[1] = radius * std::sin(theta);
        pC->assign_position( position );

        pC->custom_data["axis_a"] = semimajor;
        pC->custom_data["axis_b"] = b_axis_calc;
        pC->custom_data["axis_c"] = c_axis_calc;
        double new_volume=custom_volume_update(pC->custom_data["axis_a"], pC->custom_data["axis_b"], pC->custom_data["axis_c"]);
        pC->set_total_volume(new_volume);
        // std::cout << pC->ID << ", " <<
        theta += theta_del;
	}

    double x0 = -0.2 * Xrange;
    double xdel = 3.0;
    position[1] = -20.;
    int nbricks = 20;
    Cell_Definition* pCD = cell_definitions_by_name["wall"]; 
	for( int k=0; k <nbricks; k++ )
	{
        pC = create_cell( *pCD ); 
        pC->is_movable = false;
        position[0] = x0 + k*xdel;
        pC->assign_position( position );
	}
    position[1] = 20.;
	for( int k=0; k <nbricks; k++ )
	{
        pC = create_cell( *pCD ); 
        pC->is_movable = false;
        position[0] = x0 + k*xdel;
        pC->assign_position( position );
	}
	std::cout << std::endl; 
	
	// load cells from your CSV file (if enabled)
	load_cells_from_pugixml(); 	
	
	return; 
}

// void setup_tissue_ring( void )
// {
// 	double Xmin = microenvironment.mesh.bounding_box[0]; 
// 	double Ymin = microenvironment.mesh.bounding_box[1]; 
// 	double Zmin = microenvironment.mesh.bounding_box[2]; 

// 	double Xmax = microenvironment.mesh.bounding_box[3]; 
// 	double Ymax = microenvironment.mesh.bounding_box[4]; 
// 	double Zmax = microenvironment.mesh.bounding_box[5]; 
	
// 	if( default_microenvironment_options.simulate_2D == true )
// 	{
// 		Zmin = 0.0; 
// 		Zmax = 0.0; 
// 	}
	
// 	double Xrange = Xmax - Xmin; 
// 	double Yrange = Ymax - Ymin; 
// 	double Zrange = Zmax - Zmin; 

//     double radius = 0.4 * Xrange;
// 	Cell* pC;
	
//     std::vector<double> position = {0,0,0}; 
//     static double pi = 3.141592653589793238462643383279502884;
//     double theta_del = 2.0*pi / parameters.ints("number_of_cells");
//     // double theta = 0.0;
//     double theta = 0.2;
//     position[2] = 0.0;
// 	for( int k=0; k <parameters.ints("number_of_cells") ; k++ )
// 	{
// 		// Cell_Definition* pCD = cell_definitions_by_index[k]; 
// 		Cell_Definition* pCD = cell_definitions_by_name["ellipsey"]; 

// 		double semimajor = parameters.doubles("major_axis_2a")/2;
// 		double ecc = parameters.doubles("eccentricity");
// 		double vol = parameters.doubles("starting_cell_volume");
// 		// double b_axis_calc = pow( pow(semimajor,2)*(1-pow(ecc,2)), 0.5);
// 		// double c_axis_calc = (3*vol)/( 4*pi*pow(((1-ecc)*pow(semimajor,2) ),0.5));
//         // Didi update 7-29-22
//         double b_axis_calc = semimajor*pow( (1.0-pow(ecc,2)), 0.5);
//         double c_axis_calc = (3.0*vol)/( 4.0*pi*pow(semimajor,2)*pow((1-ecc),0.5));

// 		std::cout << "ecc " << ecc << " ... " << std::endl;
// 		std::cout << "bax " << b_axis_calc << " ... " << std::endl; 

//         // create 4 cells around the origin (in Z=0 plane)
//         pC = create_cell( *pCD ); 
//         position[0] = radius * std::cos(theta);
//         position[1] = radius * std::sin(theta);
//         pC->assign_position( position );

//         pC->custom_data["axis_a"] = semimajor;
//         pC->custom_data["axis_b"] = b_axis_calc;
//         pC->custom_data["axis_c"] = c_axis_calc;
//         double new_volume=custom_volume_update(pC->custom_data["axis_a"], pC->custom_data["axis_b"], pC->custom_data["axis_c"]);
//         pC->set_total_volume(new_volume);
//         // std::cout << pC->ID << ", " <<
//         theta += theta_del;
// 	}
// 	std::cout << std::endl; 
	
// 	// load cells from your CSV file (if enabled)
// 	load_cells_from_pugixml(); 	
	
// 	return; 
// }

std::vector<std::string> my_coloring_function( Cell* pCell )
{ return paint_by_number_cell_coloring(pCell); }


// Issy Cowlishaw
void assign_orientation(Cell* pCell, Phenotype& phenotype, double dt_)
{
    //cell.state.orientation.resize(3,0.0);
    //phenotype.state.orientation.resize(3,0.0);
    pCell->state.orientation.resize(3,0.0);
    // phenotype.state.orientation.resize(3,0.0);
    if( pCell->functions.set_orientation != NULL )
    {
        pCell->functions.set_orientation(pCell, phenotype, 0.0 );
    }
    else
    {
        //assign a random unit vector
        double theta= UniformRandom()*6.28318530717959; //rand*2*pi
        double z= 2* UniformRandom()-1;
        double temp= sqrt(1-z*z);
        pCell->state.orientation[0]= temp * cos(theta);
        pCell->state.orientation[1]= temp * sin(theta);
        pCell->state.orientation[2]= z;
    }
    return;
}
// Issy
void update_axis( Cell* pCell, Phenotype& phenotype, double dt )
{
    static double four_thirds_pi =  4.188790204786391;
    double new_vol = phenotype.volume.total;
    double scale_fac = new_vol - four_thirds_pi*pCell->custom_data["axis_a"]*pCell->custom_data["axis_b"]*pCell->custom_data["axis_c"];
    scale_fac = pow( scale_fac , 0.333333333333333333333333333333333333333 );
    pCell->custom_data["axis_a"] *= scale_fac;
    pCell->custom_data["axis_b"] *= scale_fac;
    pCell->custom_data["axis_c"] *= scale_fac;
    return;
}
// Issy
void update_motility_vector(Cell* pCell, Phenotype& phenotype, double dt_)
{
    std::cout << "------------  (custom) update_motility_vector \n";
    int a,b,angle;
    double x,y;
    double pi = 3.141592653589793238462643383279502884;
    a = pCell->custom_data["axis_a"];
    b = pCell->custom_data["axis_b"];
    angle = ((pCell->custom_data["ori_z"])*pi)/180;
    double angle90 = (90*pi)/180;
    double speed = phenotype.motility.migration_speed;
    double max_axis = std::max(a,b);
    std::cout << "max= " << max_axis;
    if (max_axis == pCell->custom_data["axis_a"])
    {
        x = speed*cos(angle);
        y = speed*sin(angle);
        phenotype.motility.motility_vector.assign(x, y);
    }
    else if (max_axis == pCell->custom_data["axis_b"])
    {
        x = speed*cos(angle+angle90);
        y = speed*sin(angle+angle90);
        phenotype.motility.motility_vector.assign(x, y);
    }
}

void phenotype_function( Cell* pCell, Phenotype& phenotype, double dt )
{ 
	// pCell->phenotype.geometry.axis_a=pCell->custom_data["axis_a"];
	// pCell->phenotype.geometry.axis_b=pCell->custom_data["axis_b"];
	// pCell->phenotype.geometry.axis_c=pCell->custom_data["axis_c"];

    if (pCell->type_name.compare("ellipsey") == 0)
    {
        pCell->custom_data["axis_a"] += 4;
        pCell->custom_data["axis_b"] += 10;
        pCell->custom_data["axis_c"] += 2;
    }
    else
    {
        double axis_length = 4.0;
        pCell->custom_data["axis_a"] = axis_length;
        pCell->custom_data["axis_b"] = axis_length;
        pCell->custom_data["axis_c"] = axis_length;
    }

    // update_motility_vector(pCell, phenotype, dt);
	
	return; 
}

void custom_function( Cell* pCell, Phenotype& phenotype , double dt )
{ return; } 

void contact_function( Cell* pMe, Phenotype& phenoMe , Cell* pOther, Phenotype& phenoOther , double dt )
{ return; } 


