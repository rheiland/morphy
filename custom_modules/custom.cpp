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
	
	// draw 1 cell 
	
	Cell* pC;
	
	for( int k=0; k < 1 ; k++ )
	{
        std::vector<double> position = {0,0,0}; 
		Cell_Definition* pCD = cell_definitions_by_index[k]; 

        // read the 

        position[2] = 0.0;

        double x0 = 70.;
        double y0 = 40.;

        // create 4 cells around the origin (in Z=0 plane)
        pC = create_cell( *pCD ); 
        position[0] = x0;
        position[1] = 0.; 
        pC->assign_position( position );
        pC->custom_data["axis_a"] = parameters.doubles("axis_a");
        pC->custom_data["axis_b"] = parameters.doubles("axis_b");
        pC->custom_data["axis_c"] = parameters.doubles("axis_c");
        double new_volume=custom_volume_update(pC->custom_data["axis_a"], pC->custom_data["axis_b"], pC->custom_data["axis_c"]);
        pC->set_total_volume(new_volume);

        pC = create_cell( *pCD ); 
        position[0] = -x0;
        position[1] = 0; 
        pC->assign_position( position );
        pC->custom_data["axis_a"] = parameters.doubles("axis_a");
        pC->custom_data["axis_b"] = parameters.doubles("axis_b");
        pC->custom_data["axis_c"] = parameters.doubles("axis_c");
        new_volume=custom_volume_update(pC->custom_data["axis_a"], pC->custom_data["axis_b"], pC->custom_data["axis_c"]);
        pC->set_total_volume(new_volume);

        pC = create_cell( *pCD ); 
        position[0] = 0;
        position[1] = y0; 
        pC->assign_position( position );
        pC->custom_data["axis_a"] = parameters.doubles("axis_a");
        pC->custom_data["axis_b"] = parameters.doubles("axis_b");
        pC->custom_data["axis_c"] = parameters.doubles("axis_c");
        new_volume=custom_volume_update(pC->custom_data["axis_a"], pC->custom_data["axis_b"], pC->custom_data["axis_c"]);
        pC->set_total_volume(new_volume);

        pC = create_cell( *pCD ); 
        position[0] = 0;
        position[1] = -y0; 
        pC->assign_position( position );
        pC->custom_data["axis_a"] = parameters.doubles("axis_a");
        pC->custom_data["axis_b"] = parameters.doubles("axis_b");
        pC->custom_data["axis_c"] = parameters.doubles("axis_c");
        new_volume=custom_volume_update(pC->custom_data["axis_a"], pC->custom_data["axis_b"], pC->custom_data["axis_c"]);
        pC->set_total_volume(new_volume);
        std::cout << "new_volume= " << new_volume << std::endl;

        // std::cout << pC->ID << ", " <<

        // std::cout << "Placing cells of type " << pCD->name << " ... " << std::endl; 
		// for( int n = 0 ; n < 1; n++ )
		// {
		// 	std::vector<double> position = {0,0,0}; 
		// 	position[0] = parameters.doubles("cx");
		// 	position[1] = parameters.doubles("cy"); 
		// 	position[2] = parameters.doubles("cz");

		// 	pC = create_cell( *pCD ); 
		// 	pC->assign_position( position );

		// 	pC->custom_data["axis_a"]=parameters.doubles("axis_a");
		// 	pC->custom_data["axis_b"]=parameters.doubles("axis_b");
		// 	pC->custom_data["axis_c"]=parameters.doubles("axis_c");
        //     std::cout << pC->ID << ", " <<
		// 	double new_volume=custom_volume_update(pC->custom_data["axis_a"], pC->custom_data["axis_b"], pC->custom_data["axis_c"]);
		// 	pC->set_total_volume(new_volume);
		// 	//std::cout<< pC->custom_data["test_list"]<<std::endl;;
		// }
	}
	std::cout << std::endl; 
	
	// load cells from your CSV file (if enabled)
	load_cells_from_pugixml(); 	
	
	return; 
}

std::vector<std::string> my_coloring_function( Cell* pCell )
{ return paint_by_number_cell_coloring(pCell); }

void phenotype_function( Cell* pCell, Phenotype& phenotype, double dt )
{ 
	pCell->phenotype.geometry.axis_a=pCell->custom_data["axis_a"];
	pCell->phenotype.geometry.axis_b=pCell->custom_data["axis_b"];
	pCell->phenotype.geometry.axis_c=pCell->custom_data["axis_c"];
	
	
	return; 
}

void custom_function( Cell* pCell, Phenotype& phenotype , double dt )
{ return; } 

void contact_function( Cell* pMe, Phenotype& phenoMe , Cell* pOther, Phenotype& phenoOther , double dt )
{ return; } 


