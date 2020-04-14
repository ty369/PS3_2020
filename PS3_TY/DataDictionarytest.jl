# ----------------------------------------------------------------------------------- #
# Copyright (c) 2018 Varnerlab
# Robert Frederick Smith School of Chemical and Biomolecular Engineering
# Cornell University, Ithaca NY 14850
# Adopted the Backbone of the code from Varnerlab
# ----------------------------------------------------------------------------------- #

function reactions(time_start,time_stop,time_step)

	# load the original data dictionary -
	data_dictionary = maximize_urea_production(time_start,time_stop,time_step)
	
	S = data_dictionary["stoichiometric_matrix"]
	S=CSV.read("/Users/ty369/Desktop/PS3_test/S_Matrix.CSV")
	# ----------------------------------------------------- #
	# Species in order:
	# 1. ATP
	# 2. L-Citrulline
	# 3. L-Aspartate
	# 4. AMP
	# 5. Diphosphate
	# 6. 2-Nomega-L-(arginino)succinate
	# 7. carbamoyl phosphate
	# 8. L-ornithine
	# 9. phosphate
	# 10. L-arginine
	# 11. H2O
	# 12. Urea
	# 13. Fumarate
	# 14. NADPH
	# 15. H+
	# 16. O2
	# 17. Nitric oxide
    # 18. NAPD+

	# repackage -
	data_dictionary["list_of_metabolite_symbols"] = list_of_metabolite_symbols
	data_dictionary["stoichiometric_matrix"] = S

	# return the modified data dictionary -
	return data_dictionary
end


function maximize_urea_production_open(time_start,time_stop,time_step)

	# load the original data dictionary -
	data_dictionary = maximize_urea_production(time_start,time_stop,time_step)

	# lets open up the side products of ec:1.14.13.39 -

	# 2: Update the reaction bounds array -
	default_flux_bounds_array = data_dictionary["default_flux_bounds_array"]

	# Vmax [mmol/gdw-hr] 15	[] --> NADPH
	default_flux_bounds_array[15,1] = 0
	default_flux_bounds_array[15,2] = 10
	# Vmax [mmol/gdw-hr] 16	[] --> H+
	default_flux_bounds_array[16,1] = 0
	default_flux_bounds_array[16,2] = 10
	# Vmax [mmol/gdw-hr] 17	[] --> Oxygen
	default_flux_bounds_array[17,1] = 0
	default_flux_bounds_array[17,2] = 10
	# Vmax [mmol/gdw-hr] 18	Nitric_oxide-->[]
	default_flux_bounds_array[18,1] = 0
	default_flux_bounds_array[18,2] = 10
	# Vmax [mmol/gdw-hr] 19	NADP --> []
	default_flux_bounds_array[19,1] = 0
	default_flux_bounds_array[19,2] = 10
	# Vmax [mmol/gdw-hr] 20	H2O --> []
	default_flux_bounds_array[20,1] = 0
	default_flux_bounds_array[20,2] = 10
	# Vmax [mmol/gdw-hr] 21	[] --> H2O
	default_flux_bounds_array[21,1] = 0
	default_flux_bounds_array[21,2] = 10

	# repackage -
	data_dictionary["default_flux_bounds_array"] = default_flux_bounds_array

	# return the updated dictionary -
	return data_dictionary
end

function maximize_urea_production(time_start,time_stop,time_step)

	# load the original data dictionary -
	data_dictionary = DataDictionary(time_start,time_stop,time_step)

	# 1: set the objective function -
	objective_coefficient_array = data_dictionary["objective_coefficient_array"]
	objective_coefficient_array[10] = -1

	# 2: Update the reaction bounds array -
	default_flux_bounds_array = data_dictionary["default_flux_bounds_array"]

	# let all exchanges go from 0,10 mmol/gDW-hr
	range_of_exchange_reactions = collect(7:21)
	for reaction_index in range_of_exchange_reactions
		default_flux_bounds_array[reaction_index,1] = 0.0
		default_flux_bounds_array[reaction_index,2] = 10.0
	end

	# don't allow water exchange -
	default_flux_bounds_array[20,2] = 0
	default_flux_bounds_array[21,2] = 0

	# we have some specific values for v1 -> v5, (values are given in the problem statment)
	E = (0.01)*(1/1000)	# mmol/gDW
	metabolic_vmax_array = [
		203*3600*E	;	# v1 ec:6.3.4.5 mmol/gDW-hr
		34.5*3600*E	;	# v2 ec:4.3.2.1 mmol/gDW-hr
		249*3600*E	;	# v3 ec:3.5.3.1 mmol/gDW-hr
		88.1*3600*E	;	# v4 ec:2.1.3.3 mmol/gDW-hr
		13.7*3600*E	;	# v5 ec:1.14.13.39 mmol/gDW-hr
		13.7*3600*E	;	# v6 ec:1.14.13.39 mmol/gDW-hr
	]
	

	range_of_cycle_reactions = collect(1:6)
	for reaction_index in range_of_cycle_reactions
		default_flux_bounds_array[reaction_index,1] >= 0.0
		default_flux_bounds_array[reaction_index,2] <= metabolic_vmax_array[reaction_index]
	end

	# repackage -
	data_dictionary["default_flux_bounds_array"] = default_flux_bounds_array
	data_dictionary["objective_coefficient_array"] = objective_coefficient_array

	# return the updated dictionary -
	return data_dictionary
end


# ----------------------------------------------------------------------------------- #
# Function Originally Generated on: 2018-03-15T00:00:56.939
# Updated: 2020-04-14
#
# Input arguments:
# time_start::Float64 => Simulation start time value (scalar)
# time_stop::Float64 => Simulation stop time value (scalar)
# time_step::Float64 => Simulation time step (scalar)
#
# Output arguments:
# data_dictionary::Dict{AbstractString,Any} => Dictionary holding model and simulation parameters as key => value pairs
# ----------------------------------------------------------------------------------- #
function DataDictionary(time_start,time_stop,time_step)

	# Load the stoichiometric network from disk -
	stoichiometric_matrix =readdlm("Elements.txt")
	(number_of_species,number_of_reactions) = size(stoichiometric_matrix)
	E = (0.01)*(1/1000)	# mmol/gDW
	metabolic_vmax_array = [
	203*3600*E	;	# v1 ec:6.3.4.5 mmol/gDW-hr
	34.5*3600*E	;	# v2 ec:4.3.2.1 mmol/gDW-hr
	249*3600*E	;	# v3 ec:3.5.3.1 mmol/gDW-hr
	88.1*3600*E	;	# v4 ec:2.1.3.3 mmol/gDW-hr
	13.7*3600*E	;	# v5 ec:1.14.13.39 mmol/gDW-hr
	13.7*3600*E	;	# v6 ec:1.14.13.39 mmol/gDW-hr
	]
	# Modify vmax based on km and [S] from Park et al's paper, v=kcat*[E]*([S]/([S]+[Km]))
	v1=metabolic_vmax_array[1]*0.92256*0.98977
	v2=metabolic_vmax_array[2]*1
	v3=metabolic_vmax_array[3]*0.1421*1
	v4=metabolic_vmax_array[4]*0.2627
	v5=metabolic_vmax_array[5]*0.98646*1
	v5r=metabolic_vmax_array[6]*(1)
	
	
	# Setup default flux bounds array -
	default_bounds_array = [
		0	v1	;	# Vmax [mmol/gdw-hr] 1	ATP+L-Citrulline+L-Aspartate-->AMP+Diphosphate+N-(L-Arginino)succinate
		0	v2;	# Vmax [mmol/gdw-hr] 2	N-(L-Arginino)succinate-->Fumarate+L-Arginine
		0	v3;	# Vmax [mmol/gdw-hr] 3	L-Arginine+H2O-->L-Ornithine+Urea
		0	v4;	# Vmax [mmol/gdw-hr] 4	Carbamoyl_phosphate+L-Ornithine-->Orthophosphate+L-Citrulline
		0	v5;	# Vmax [mmol/gdw-hr] 5	2.0*L-Arginine+4.0*Oxygen+3.0*NADPH+3.0*H-->2.0*Nitric_oxide+2.0*L-Citrulline+3.0*NADP+4.0*H2O
		0	v5r;	# Vmax [mmol/gdw-hr] 6	2.0*Nitric_oxide+2.0*L-Citrulline+3.0*NADP+4.0*H2O-->2.0*L-Arginin+4.0*Oxygen+3.0*NADPH+3.0*H
		0	2.2148976	;	# Vmax [mmol/gdw-hr] 7	[] -->Carbamoyl_phosphate
		0	2.2148976	;	# Vmax [mmol/gdw-hr] 8	[] -->L-Aspartate
		0	2.2148976	;	# Vmax [mmol/gdw-hr] 9	Fumarate--> []
		0	2.2148976	;	# Vmax [mmol/gdw-hr] 10	Urea--> []
		0	2.2148976	;	# Vmax [mmol/gdw-hr] 11	[] --> ATP
		0	2.2148976	;	# Vmax [mmol/gdw-hr] 12	AMP --> []
		0	2.2148976	;	# Vmax [mmol/gdw-hr] 13	Diphosphate--> []
		0	2.2148976	;	# Vmax [mmol/gdw-hr] 14	Orthophosphate--> []
		0	2.2148976	;	# Vmax [mmol/gdw-hr] 17	[] -->Oxygen
		0	2.2148976	;	# Vmax [mmol/gdw-hr] 15	[] -->NADPH
		0	2.2148976	;	# Vmax [mmol/gdw-hr] 16	[] --> H
		0	2.2148976	;	# Vmax [mmol/gdw-hr] 18	Nitric_oxide--> []
		0	2.2148976	;	# Vmax [mmol/gdw-hr] 19	NADP --> []
		0	2.2148976	;	# Vmax [mmol/gdw-hr] 20	H2O--> []
		0	2.2148976	;	# Vmax [mmol/gdw-hr] 21	[] --> H2O
	];


	# Setup default species bounds array -
	species_bounds_array = [
		0.0	0.0	;	# 1. ATP
		0.0	0.0	;	# 2. L-Citrulline
		0.0	0.0	;	# 3. L-Aspartate
		0.0	0.0	;	# 4. AMP
		0.0	0.0	;	# 5. Diphosphate
		0.0	0.0	;	# 6. 2-Nomega-L-(arginino)succinate
		0.0	0.0	;	# 7. carbamoyl phosphate
		0.0	0.0	;	# 8. L-ornithine
		0.0	0.0	;	# 9. phosphate
		0.0	0.0	;	# 10. L-arginine
		0.0	0.0	;	# 11. H2O
		0.0	0.0	;	# 12. Urea
		0.0	0.0	;	# 13. Fumarate
		0.0	0.0	;	# 14. NADPH
		0.0	0.0	;	# 15. H+
		0.0	0.0	;	# 16. O2
		0.0	0.0	;	# 17. Nitric oxide
		0.0	0.0	;	# 18. NAPD+
	];

	# Min/Max flag - default is minimum -
	is_minimum_flag = true

	# Setup the objective coefficient array - v1-v5, b1-b15 corespond to the ones in excel
	objective_coefficient_array = [

		0.0	;	# 1 v1
		0.0	;	# 2 v2
		0.0	;	# 3 v3
		0.0	;	# 4 v4
		0.0	;	# 5 v5
		0.0	;	# 6 v5_reverse
		0.0	;	# 7 b1
		0.0	;	# 8 b2
		0.0	;	# 9 b3
		0.0	;	# 10 b4
		0.0	;	# 11 b5
		0.0	;	# 12 b6
		0.0	;	# 13 b7
		0.0	;	# 14 b8
		0.0	;	# 15 b9
		0.0	;	# 16 b10
		0.0	;	# 17 b11
		0.0	;	# 18 b12
		0.0	;	# 19 b13
		0.0	;	# 20 b14
		0.0	;	# 21 b14_reverse:

	];

	# List of reation strings - used to write flux report
	list_of_reaction_strings = [
		"v1::ATP+L-Citrulline+L-Aspartate-->AMP+Diphosphate+N-(L-Arginino)succinate"	;	# 1
		"v2::N-(L-Arginino)succinate-->Fumarate+L-Arginine"	;	# 2
		"v3::L-Arginine+H2O-->L-Ornithine+Urea"	;	# 3
		"v4::Carbamoyl_phosphate+L-Ornithine-->Orthophosphate+L-Citrulline"	;	# 4
		"v5::2.0*L-Arginine+4.0*Oxygen+3.0*NADPH+3.0*H-->2.0*Nitric_oxide+2.0*L-Citrulline+3.0*NADP+4.0*H2O"	;	# 5
		"v5_reverse::2.0*Nitric_oxide+2.0*L-Citrulline+3.0*NADP+4.0*H2O-->2.0*L-Arginine+4.0*Oxygen+3.0*NADPH+3.0*H"	;	# 6
		"b1::[]-->Carbamoyl_phosphate"	;	# 7
		"b2::[]-->L-Aspartate"	;	# 8
		"b3::Fumarate-->[]"	;	# 9
		"b4::Urea-->[]"	;	# 10
		"b5::[]-->ATP"	;	# 11
		"b6::AMP-->[]"	;	# 12
		"b7::Diphosphate-->[]"	;	# 13
		"b8::Orthophosphate-->[]"	;	# 14
		"b11::[]-->Oxygen"	;	# 15
		"b9::[]-->NADPH"	;	# 16
		"b10::[]-->Hc"	;	# 17
		"b12::Nitric_oxide-->[]"	;	# 18
		"b13::NADP-->[]"	;	# 19
		"b14::H2O-->[]"	;	# 20
		"b14_reverse::[]-->H2O"	;	# 21
	];

	# List of metabolite strings - used to write flux report
	list_of_metabolite_symbols = [
		"AMP"	;	# 1
		"ATP"	;	# 2
		"Carbamoyl_phosphate"	;	# 3
		"Diphosphate"	;	# 4
		"Fumarate"	;	# 5
		"H2O"	;	# 6
		"H"	;	# 7
		"L-Arginine"	;	# 8
		"L-Aspartate"	;	# 9
		"L-Citrulline"	;	# 10
		"L-Ornithine"	;	# 11
		"N-(L-Arginino)succinate"	;	# 12
		"NADPH"	;	# 13
		"NADP"	;	# 14
		"Nitric_oxide"	;	# 15
		"Orthophosphate"	;	# 16
		"Oxygen"	;	# 17
		"Urea"	;	# 18
	];


	# Metabolic Vmax array (units: mmol/B-hr) -
	metabolic_vmax_array = [
		v1	;	# Vmax [mmol/gdw-hr] 1	v1
		v2	;	# Vmax [mmol/gdw-hr] 2	v2
		v3	;	# Vmax [mmol/gdw-hr] 3	v3
		v4	;	# Vmax [mmol/gdw-hr] 4	v4
		v5	;	# Vmax [mmol/gdw-hr] 5	v5
		v5r	;	# Vmax [mmol/gdw-hr] 6	v5r
		2.2148976	;	# Vmax [mmol/gdw-hr] 7	b1
		2.2148976	;	# Vmax [mmol/gdw-hr] 8	b2
		2.2148976	;	# Vmax [mmol/gdw-hr] 9	b3
		2.2148976	;	# Vmax [mmol/gdw-hr] 10	b4
		2.2148976	;	# Vmax [mmol/gdw-hr] 11	b5
		2.2148976	;	# Vmax [mmol/gdw-hr] 12	b6
		2.2148976	;	# Vmax [mmol/gdw-hr] 13	b7
		2.2148976	;	# Vmax [mmol/gdw-hr] 14	b8
		2.2148976	;	# Vmax [mmol/gdw-hr] 15	b11
		2.2148976	;	# Vmax [mmol/gdw-hr] 16	b9
		2.2148976	;	# Vmax [mmol/gdw-hr] 17	b10
		2.2148976	;	# Vmax [mmol/gdw-hr] 18	b12
		2.2148976	;	# Vmax [mmol/gdw-hr] 19	b13
		2.2148976	;	# Vmax [mmol/gdw-hr] 20	b20
		2.2148976	;	# Vmax [mmol/gdw-hr] 21	b21
	];

	# Metabolic saturation constant array (units mM) -
	number_of_metabolic_rates = length(metabolic_vmax_array)
	metabolic_saturation_constant_array = 0.130*ones(number_of_metabolic_rates*number_of_species)

	# ------------------------------------------------------------------------------------------#
	# constants (from bionumbers)       units
	# ------------------------------------------------------------------------------------------#
	cell_diameter = 1.1                             # mum
	mass_of_single_cell = 2.8e-13                   # g
	number_of_rnapII = 4600            	            # copies/cells
	number_of_ribosome = 50000         	            # copies/cells
	mRNA_half_life_TF = 0.083                       # hrs
	protein_half_life = 70                          # hrs
	infrastructure_half_life = 300					# hrs
	doubling_time_cell = 0.33                       # hrs
	max_translation_rate = 16.5                     # aa/sec
	max_transcription_rate = 60.0                   # nt/sec
	transcription_initiation_time_contstant = 400  # sec
	average_transcript_length = 1200   	            # nt
	average_protein_length = 400       	            # aa
	fraction_nucleus = 0.0             	            # dimensionless
	av_number = 6.02e23                             # number/mol
	avg_gene_number = 2                             # number of copies of a gene
	polysome_number = 4					            # number of ribsomoses per transcript
	# ------------------------------------------------------------------------------------------#
	#
	# ------------------------------------------------------------------------------------------#
	# Calculate constants using bionumber values
	# ------------------------------------------------------------------------------------------#
	# Calculate the volume (convert to L)
	V = ((1-fraction_nucleus)*(1/6)*(3.14159)*(cell_diameter)^3)*(1e-15)

	# Calculate the rnapII_concentration and ribosome_concentration
	rnapII_concentration = number_of_rnapII*(1/av_number)*(1/mass_of_single_cell)*1e6       # mumol/gdw
	ribosome_concentration = number_of_ribosome*(1/av_number)*(1/mass_of_single_cell)*1e6   # mumol/ddw

	# degrdation rate constants -
	degradation_constant_mRNA = -(1/mRNA_half_life_TF)*log(0.5)                           # hr^-1
	degradation_constant_protein = -(1/protein_half_life)*log(0.5)                        # hr^-1
	degrdation_constant_infrastructure = -(1/infrastructure_half_life)*log(0.5)			# hr^-1

	# kcats for transcription and translation -
	kcat_transcription = max_transcription_rate*(3600/average_transcript_length)            # hr^-1
	kcat_translation = polysome_number*max_translation_rate*(3600/average_protein_length)   # hr^-1

	# kcat for transcription initiation -
	kcat_transcription_initiation = ((1/3600)*transcription_initiation_time_contstant)^-1   # hr^-1
	kcat_translation_initiation = 10*kcat_transcription_initiation                          # hr^-1

	# Maximum specific growth rate -
	maximum_specific_growth_rate = (1/doubling_time_cell)*log(2)                          # hr^-1

	# What is the average gene concentration -
	avg_gene_concentration = avg_gene_number*(1/mass_of_single_cell)*(1/V)*1e9              # nmol/gdw

	# How fast do my cells die?
	death_rate_constant = 0.05*maximum_specific_growth_rate                                 # hr^-1

	# Saturation constants for translation and trascription -
	saturation_transcription = 4600*(1/av_number)*(1/mass_of_single_cell)*1e9               # nmol/gdw
	saturation_translation = 150000*(1/av_number)*(1/mass_of_single_cell)*1e6               # nmol/gdw
	# -------------------------------------------------------------------------------------------#


	# Alias the txtl parameters -
	txtl_parameter_dictionary = Dict{AbstractString,Any}()
	txtl_parameter_dictionary["rnapII_concentration"] = rnapII_concentration  # muM
	txtl_parameter_dictionary["ribosome_concentration"] = ribosome_concentration # muM
	txtl_parameter_dictionary["degradation_constant_mRNA"] = degradation_constant_mRNA  # hr^-1
	txtl_parameter_dictionary["degradation_constant_protein"] = degradation_constant_protein  # hr^-1
	txtl_parameter_dictionary["degrdation_constant_infrastructure"] = degrdation_constant_infrastructure  # hr^-1
	txtl_parameter_dictionary["kcat_transcription"] = kcat_transcription  # hr^-1
	txtl_parameter_dictionary["kcat_translation"] = kcat_translation  # hr^-1
	txtl_parameter_dictionary["maximum_specific_growth_rate"] = maximum_specific_growth_rate  # hr^-1
	txtl_parameter_dictionary["death_rate_constant"] = death_rate_constant
	txtl_parameter_dictionary["avg_gene_concentration"] = avg_gene_concentration
	txtl_parameter_dictionary["saturation_constant_transcription"] = saturation_transcription
	txtl_parameter_dictionary["saturation_constant_translation"] = saturation_translation
	txtl_parameter_dictionary["average_transcript_length"] = average_transcript_length
	txtl_parameter_dictionary["average_protein_length"] = average_protein_length
	txtl_parameter_dictionary["kcat_transcription_initiation"] = kcat_transcription_initiation
	txtl_parameter_dictionary["kcat_translation_initiation"] = kcat_translation_initiation

	# Setup species abundance array -
	species_abundance_array = [
		0.0	;	# 1 AMP
		0.0	;	# 2 ATP
		0.0	;	# 3 Carbamoyl_phosphate
		0.0	;	# 4 Diphosphate
		0.0	;	# 5 Fumarate
		0.0	;	# 6 H2O
		0.0	;	# 7 H
		0.0	;	# 8 L-Arginine
		0.0	;	# 9 L-Aspartate
		0.0	;	# 10 L-Citrulline
		0.0	;	# 11 L-Ornithine
		0.0	;	# 12 N-(L-Arginino)succinate
		0.0	;	# 13 NADPH
		0.0	;	# 14 NADP
		0.0	;	# 15 Nitric_oxide
		0.0	;	# 16 Orthophosphate
		0.0	;	# 17 Oxygen
		0.0	;	# 18 Urea
	];

	# =============================== DO NOT EDIT BELOW THIS LINE ============================== #
	data_dictionary = Dict{AbstractString,Any}()
	data_dictionary["stoichiometric_matrix"] = stoichiometric_matrix
	data_dictionary["objective_coefficient_array"] = objective_coefficient_array
	data_dictionary["default_flux_bounds_array"] = default_bounds_array;
	data_dictionary["species_abundance_array"] = species_abundance_array
	data_dictionary["species_bounds_array"] = species_bounds_array
	data_dictionary["list_of_reaction_strings"] = list_of_reaction_strings
	data_dictionary["list_of_metabolite_symbols"] = list_of_metabolite_symbols
	data_dictionary["is_minimum_flag"] = is_minimum_flag
	data_dictionary["number_of_species"] = number_of_species
	data_dictionary["number_of_reactions"] = number_of_reactions
	data_dictionary["metabolic_saturation_constant_array"] = metabolic_saturation_constant_array
	data_dictionary["metabolic_vmax_array"] = metabolic_vmax_array
	data_dictionary["characteristic_enzyme_abundance_mM"] = 4.7326871794871794e-5
	data_dictionary["volume_of_cell"] = V
	data_dictionary["mass_of_single_cell"] = mass_of_single_cell
	data_dictionary["txtl_parameter_dictionary"] = txtl_parameter_dictionary
	# =============================== DO NOT EDIT ABOVE THIS LINE ============================== #
	return data_dictionary
end
