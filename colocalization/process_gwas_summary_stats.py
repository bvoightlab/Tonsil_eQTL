'''
Testing Variables

config_file="C:/Users/mitch/Documents/UPenn/Voight_Lab/Trait_Colocalization/config_files/summary_stats_cleaning.testing.yml"

trait_file_loc="C:/Users/mitch/Documents/UPenn/Voight_Lab/Trait_Colocalization/Deprecated/test.txt.gz"

conversion_dict = {'GRCh37_to_GRCh38': "C:/Users/mitch/Documents/UPenn/Voight_Lab/Trait_Colocalization/Deprecated/GRCh37_to_GRCh38.json", 
                   'GRCh37_to_rsid': "C:/Users/mitch/Documents/UPenn/Voight_Lab/Trait_Colocalization/Deprecated/GRCh37_to_rsid.json", 
                   'GRCh38_to_GRCh37': "C:/Users/mitch/Documents/UPenn/Voight_Lab/Trait_Colocalization/Deprecated/GRCh38_to_GRCh37.json", 
                   'GRCh38_to_rsid': "C:/Users/mitch/Documents/UPenn/Voight_Lab/Trait_Colocalization/Deprecated/GRCh38_to_rsid.json",
                   'rsid_to_GRCh37': "C:/Users/mitch/Documents/UPenn/Voight_Lab/Trait_Colocalization/Deprecated/rsid_to_GRCh37.json", 
                   'rsid_to_GRCh38': "C:/Users/mitch/Documents/UPenn/Voight_Lab/Trait_Colocalization/Deprecated/rsid_to_GRCh38.json"}
dbsnp_vcf_loc_dict =  {'GRCh37': "C:/Users/mitch/Documents/UPenn/Voight_Lab/Trait_Colocalization/Deprecated/dbSNP155.GRCh37.abridged.vcf.gz", 
                   'GRCh38': "C:/Users/mitch/Documents/UPenn/Voight_Lab/Trait_Colocalization/Deprecated/dbSNP155.GRCh38.abridged.vcf.gz"}


with open("C:/Users/mitch/Documents/UPenn/Voight_Lab/Trait_Colocalization/Deprecated/col_out_size_dict_example.json") as json_file:
            col_out_size_dict = json.load(json_file)
final_alignment='GRCh38'
min_pval=1e-300
all_aligns=['GRCh37', 'GRCh38']
retain_align=True
conversion_dict = {'GRCh37_to_GRCh38': "C:/Users/mitch/Documents/UPenn/Voight_Lab/Trait_Colocalization/Deprecated/GRCh37_to_GRCh38.json", 
                   'GRCh37_to_rsid': "C:/Users/mitch/Documents/UPenn/Voight_Lab/Trait_Colocalization/Deprecated/GRCh37_to_rsid.json", 
                   'GRCh38_to_GRCh37': "C:/Users/mitch/Documents/UPenn/Voight_Lab/Trait_Colocalization/Deprecated/GRCh38_to_GRCh37.json", 
                   'GRCh38_to_rsid': "C:/Users/mitch/Documents/UPenn/Voight_Lab/Trait_Colocalization/Deprecated/GRCh38_to_rsid.json",
                   'rsid_to_GRCh37': "C:/Users/mitch/Documents/UPenn/Voight_Lab/Trait_Colocalization/Deprecated/rsid_to_GRCh37.json", 
                   'rsid_to_GRCh38': "C:/Users/mitch/Documents/UPenn/Voight_Lab/Trait_Colocalization/Deprecated/rsid_to_GRCh38.json"}

'''

#import modules that will be needed
import sys #Taking command line parameters as inputs
import gzip #Opening gzipped input files
import json #For writing the vcf dictionaries to file for future use
import yaml
from os import path
import getopt

#Define help function to be called if there is a problem
def help(exit_num=1):
    print("""-----------------------------------------------------------------
ARGUMENTS
    -c => <yaml> config file holding all necessary information REQUIRED
ASSUMPTIONS
    * all traits to be tested must be listed in the .yaml configuration file
    * the column headers of each data type must be listed for each trait
    * each trait must be listed as case/control (cc) or continuous (quant)
    * missing data should be coded with an NA
    * if traits are missing critical data dictionaries must be input
    * if dictionaries are not yet created, a dbSNP vcf must be provided
    * Include column separator in config file if files are not delimted by ",", "\t", " ", or "|"
""")
    sys.exit(exit_num)

#Define main function that assesses the input flags and returns help message if needed
def main(argv): 
    try: 
        opts, args = getopt.getopt(sys.argv[1:], "c:nh")
                                                              
    except getopt.GetoptError:
        print("ERROR: Incorrect usage of getopts flags!")
        help()

    options_dict = dict(opts)
    
    ## Check for help flag and return empty string to break function
    help_ind = options_dict.get('-h', False)
    if help_ind != False:
        help()
    
    ## Required arguments
    try:
        config_file = options_dict['-c']
    except KeyError:
        print("ERROR: Incorrect usage of getopts flags!")
        help()
        
    #Make a one-time call to the driver script
    driver(config_file)        

#This is the driver function that will drive the script forward and is called by
#the main function after the dictionary has been read in
def driver(config_file):
    
    #Denote list of key words to denote missing info
    miss_denoters = ['NA', 'None', "Null", "na", "N/a", "Na", "n/a", "Missing", "Miss", "MISS", "missing", "miss", "none"]
    #Create an error_checker to see if any error messages popped up
    #This will allow multiple error messages to populate during a failed run allowing multiple fixes
    error_checker = False
    #Load the yaml configuration dictionary and extract subdicts/locations
    try:
        config_dict = yaml.load(open(config_file, 'r'), Loader=yaml.SafeLoader)
    except (yaml.scanner.ScannerError, yaml.parser.ParserError):
        print("ERROR: Config File is not readable. Please review yaml format and template.")
        error_checker = True
        
    ###########################################################################################################
    ################# This section of code checks for obvious errors in the config file #######################
    ###########################################################################################################
    #Check if a traits dictionary was given
    try:
        traits_dict = config_dict['traits']
    except KeyError:
        print("ERROR: No traits dictionary provided or config file formatted incorrectly")
        error_checker = True
    #Check if an output_dir was provided
    try:
        output_dir = config_dict['output_dir']
    except KeyError:
        print("ERROR: No output directory provided or config file formatted incorrectly")
        error_checker = True
    #Check if a final alignment was given
    try:
        final_alignment = config_dict['final_alignment']
    except:
        print("ERROR: No output directory provided or config file formatted incorrectly")
        error_checker = True
    #Check if conversion dictionaries were given
    try: 
        conversion_dict = config_dict['conversion_dictionary_locs']
    except KeyError:
        conversion_dict = {}
    #Check if dbsnp_vcf_locs were given
    try: 
        dbsnp_vcf_loc_dict = config_dict['dbsnp_vcf_file_locs']
    except KeyError:
        dbsnp_vcf_loc_dict = {}
    #Check if retain_all_alignments flag is present
    try: 
        retain_align = config_dict['retain_all_alignments']
        #Correct to boolean if boolean not detected
        if type(retain_align) != bool:
            if type(retain_align) == str:
                if retain_align.lower() == 'true':
                    retain_align = True
                elif retain_align.lower() == 'false':
                    retain_align = False
                else:
                    error_checker = True
                    print("ERROR: Unrecognized option selected for retain-all-alignments flag, select True/False")
            elif type(retain_align) == int:
                if retain_align == 1:
                    retain_align = True
                elif retain_align == 0:
                    retain_align = False
                else: 
                    error_checker = True
                    print("ERROR: Unrecognized option selected for retain-all-alignments flag, select True/False")
            else: 
                error_checker = True
                print("ERROR: Unrecognized option selected for retain-all-alignments flag, select True/False")
    except KeyError:
        retain_align = False
    #Break at this point if error_checker == True
    if error_checker == True:
        sys.exit(1)
    #Check the alignments to see which alignments are specified and also check whether rsid/position included
    align_rsid_pos = {}
    for trait in traits_dict.keys():
        #Check alignment
        try:
            align_rsid_pos[trait] = [traits_dict[trait]["align"]]
        except KeyError:
            align_rsid_pos[trait] = ['NA']
        #Add position and rsid info onto the align_rsid_pos value
        try:
            align_rsid_pos[trait].extend([traits_dict[trait]["rsid"], traits_dict[trait]["pos"]])
            #Convert any missing values to None
            if align_rsid_pos[trait][0] in miss_denoters:
               align_rsid_pos[trait][0] = None
            if align_rsid_pos[trait][1] in miss_denoters:
               align_rsid_pos[trait][1] = None
            if align_rsid_pos[trait][2] in miss_denoters:
               align_rsid_pos[trait][2] = None
        except KeyError:
            print("ERROR: rsid and/or pos not specified in config file for trait, " + trait)
            error_checker = True
    #Break at this point if error_checker == True
    if error_checker == True:
        sys.exit(1)
        
    ###########################################################################################################
    ################# This section of code checks for errors due to missing dictionaries ######################
    ###########################################################################################################
    
    #Make a dictionary called conversion_dict_need to track which provided dictionaries are necessary and only make them
    conversion_need_dict = {}
    for conversion in conversion_dict.keys():
        conversion_need_dict[conversion] = False
    #Get list of all alignments needed and those that are in input files
    all_inp_aligns = list(set([item[0] for item in list(align_rsid_pos.values())]))
    all_aligns = all_inp_aligns[:]
    all_aligns.append(final_alignment)
    all_aligns = list(set(all_aligns))
    all_inp_aligns = list(set(all_inp_aligns))
    #Remove any Nones from the all_aligns and all_inp_aligns list - this occurs when pos is missing and alignment is not specified
    if None in all_aligns:
        all_aligns.remove(None)
    if None in all_inp_aligns:
        all_inp_aligns.remove(None)
    #If retain align flag is True then
    #Check if interconversion dictionaries exist for all alignment conversions (except the self conversion)
    if retain_align == True:
        for i in range(len(all_inp_aligns)):
            for j in range(len(all_aligns)):
                if i != j:
                    if all_aligns[i] + "_to_" + all_aligns[j] not in conversion_dict.keys():
                        print("ERROR: Retain all alignments option selected but " + all_aligns[i] + "_to_" + all_aligns[j] + " dictionary is not provided")
                        error_checker = True
                    else:
                        #Update conversion need dict to reflect the dictionary need
                        conversion_need_dict[all_aligns[i] + "_to_" + all_aligns[j]] = True
    #Cycle through all the traits in the align_rsid_pos dictionary to verify that all the dictionaries are present
    #Need to do it as separate loop to previously consider impact due to flag for all alignment being flipped
    for trait in align_rsid_pos.keys():
        if align_rsid_pos[trait][2] == None:
            #Case where position is not supplied
            if align_rsid_pos[trait][1] == None:
                #Neither SNP nor position is supplied so throw error
                print("ERROR: rsid and position are both missing from input file for trait, " + trait)
                error_checker = True
            elif "rsid_to_" + final_alignment not in conversion_dict.keys():
                #rsid is missing and no dictionary is provided for conversion to final alignment
                print("ERROR: rsid_to_" + final_alignment + " dictionary is needed but not specified for trait, " + trait)
                error_checker = True
            else: 
                #Check Case where retain_align == True
                if retain_align == True:
                    for align in all_aligns:
                        #Update need dict to reflect dictionary need for all alignments
                        conversion_need_dict["rsid_to_" + align] = True
                else:
                    #Update conversion need dict to reflect the dictionary need for final alignment only
                    conversion_need_dict["rsid_to_" + final_alignment] = True
        elif align_rsid_pos[trait][0] == None: 
            #Case where alignment is not supplied
            if align_rsid_pos[trait][2] != None:
                #Alignment is not provided, but position is, so throw error
                print("ERROR: alignment is not provided for trait, " + trait + ", but position is - need to specify alignment") 
                error_checker = True
        elif align_rsid_pos[trait][1] == None:
            #Case where rsid is missing but position and alignment have been supplied
            #Check if alignment == final_alignment and verify alignment conversion dictionary if not
            if align_rsid_pos[trait][0] != final_alignment and align_rsid_pos[trait][0] + "_to_" + final_alignment not in conversion_dict.keys():
                #Dictionary is missing but needed so throw error
                print("ERROR: " + align_rsid_pos[trait][0] + "_to_" + final_alignment + " dictionary is needed but not specified for trait, " + trait)
                error_checker = True
            elif align_rsid_pos[trait][0] != final_alignment: 
                #Update need for the alignment to final_alignment dictionary
                conversion_need_dict[align_rsid_pos[trait][0] + "_to_" + final_alignment] = True
            #Since SNP is missing verify that pos to rsid dictionary is specified
            if align_rsid_pos[trait][0] + "_to_rsid"  not in conversion_dict.keys():
                #No dictionary specified for finding rsid so throw error
                print("ERROR: " + align_rsid_pos[trait][0] + "_to_rsid dictionary is needed for trait, " + trait + ", but not specified")
                error_checker = True
            else:
                #Update need for for alignment to rsid dictionary
                conversion_need_dict[align_rsid_pos[trait][0] + "_to_rsid"] = True
        else:
            #Case where alignment, rsid, and pos are all given. We only need to verify if alignment dictionaries are needed
            if align_rsid_pos[trait][0] != final_alignment and align_rsid_pos[trait][0] + "_to_" + final_alignment not in conversion_dict.keys():
                #Missing alignment conversion dictionary
                print("ERROR: " + align_rsid_pos[trait][0] + "_to_" + final_alignment + " dictionary is needed but not specified for trait, " + trait)
                error_checker = True
            elif align_rsid_pos[trait][0] != final_alignment:
                #Case where the alignment doesn't match the final alignment - we need to update dictionary need
                conversion_need_dict[align_rsid_pos[trait][0] + "_to_" + final_alignment] = True
    #Check if mininmum p-value is given otherwise use 1e-300
    try: 
        min_pval = float(config_dict['min_pval'])
    except KeyError:
        min_pval = float(1e-300)
        print("WARNING: No minimum p-value input, so defaulting to 1e-300")
        
    ######################################################################################################
    ######################################################################################################
    ######################################################################################################
    
    #Call handle_conversion for the dictionaries listed in the conversion_dict
    print("NOTE: Beginning dictionary checking/creation")
    no_conversion_probs = True
    for conversion_type in conversion_need_dict.keys(): 
        no_conversion_probs = no_conversion_probs * handle_conversion(conversion_dict, conversion_need_dict, dbsnp_vcf_loc_dict, conversion_type)
    #If no problems encountered in dictionary creation/processing than begin creating the final trait files
    #First create bed_no_errors to check for errors
    bed_no_errors = True
    if no_conversion_probs == True:
        print("NOTE: Dictionaries processed successfully. Beginning bed file creation.")
        for trait in traits_dict.keys():
            #Create dict for the worked trait
            worked_trait_dict = traits_dict[trait]
            #Call out_bed_trait function (Make multiprocess later?)
            col_out_size_dict = get_orig_col_and_sizes(trait, worked_trait_dict, align_rsid_pos, output_dir)
            #Check to see if an error occured during the col_out_size_dict creation
            if col_out_size_dict != False:
                bed_no_errors = bed_no_errors * out_bed_trait(trait, final_alignment, min_pval, retain_align, all_aligns, col_out_size_dict, conversion_dict)
            else: 
                bed_no_errors = False
        #Check if still no errors 
        if bed_no_errors == True:
            #Reaching this end point without encountering an error implies that all files were successfully created!
            print('SUCCESS: All bed files generated and output to ' + output_dir)
    else:
        print("ERROR: No files created due to issues with processing field conversions. Please review config file.")
        sys.exit(1)

#Create function that will check whether a dictionary exists, if it needs to, and make it if so
#Works for a specific dictionary type and returns a boolean reflecting if errors occured
def handle_conversion(conversion_dict, conversion_need_dict, dbsnp_vcf_loc_dict, conversion_type):
    #Check whether dictionaries are needed and 
    #Create them if needed after verifying vcf exists.
    if conversion_need_dict[conversion_type] == False:
        print('WARNING: Location listed for ' + conversion_type + ' dictionary, but no conversion needed.')
        return True
    elif conversion_dict == {}:
        print("ERROR: No conversion jsons provided, but " + conversion_type + "dictionary is needed.")
        print("If dbsnp vcf is provided, a location for the jsons to be populated into is still needed.")
        return False
    elif path.exists(conversion_dict[conversion_type]):
        print("NOTE: " + conversion_type + " conversion dictionary already exists and will be used.")
        return True
    #If we've gotten to this point it means that the conversion dict is needed and will have to be created
    #Check if we can identify a to and from part of the conversion
    if conversion_type.find("_to_") != -1:
        from_type = conversion_type[:conversion_type.find("_to_")]
        to_type = conversion_type[conversion_type.find("_to_") + 4:]
    else: 
        print("ERROR: Incorrectly formatted dictionary name. Dictionary should be of format A_to_B")
        return False
    #Check that any vcfs were supplied
    if dbsnp_vcf_loc_dict != {}:
        #Check if either end of conversion is rsid and call appropriate function
        if from_type == "rsid":
            #Attempt to denote vcf loc
            try:
                dbsnp_vcf_file_loc = dbsnp_vcf_loc_dict[to_type]
            except KeyError:
                print("ERROR: The " + to_type + " vcf was not provided and is needed")
                return False
            return_val = make_rsid_to_pos_dict(to_type, dbsnp_vcf_file_loc, conversion_dict[conversion_type])
            return return_val #Will be False if an error occured/Error message inside function
        elif to_type == "rsid":
            #Attempt to denote vcf loc
            try:
                dbsnp_vcf_file_loc = dbsnp_vcf_loc_dict[from_type]
            except KeyError:
                print("ERROR: The " + from_type + " vcf was not provided but is needed")
                return False
            return_val = make_pos_to_rsid_dict(from_type, dbsnp_vcf_file_loc, conversion_dict[conversion_type])
            return return_val #Will be False if an error occured/Error message inside function
        else:
            #Attempt to denote vcf locs
            try:
                dbsnp_vcf_from_loc = dbsnp_vcf_loc_dict[from_type]
                dbsnp_vcf_to_loc = dbsnp_vcf_loc_dict[to_type]
            except KeyError:
                print("ERROR: The " + from_type + " vcf and/or " + to_type + " vcf was not provided but both are needed")
                return False
            return_val = make_pos_to_pos_dict(from_type, to_type, dbsnp_vcf_from_loc, dbsnp_vcf_to_loc, conversion_dict[conversion_type])
            return return_val #Will be False if an error occured/Error message inside function
    else:
        print("ERROR: " + conversion_type + " dictionary does not exist and is needed but dbsnp vcf(s) not provided in config file")
        return False

#Create a function that will take in the dbsnp_vcf_file_loc and make a dictionary
#for rsid_to_pos 
def make_rsid_to_pos_dict(to_type, dbsnp_vcf_file_loc, dict_loc):
    
    #Open the input file but don't read in yet
    from itertools import islice
    try:
        dbsnp_vcf_file = opener(dbsnp_vcf_file_loc)
    except FileNotFoundError:
        print("ERROR: File loc for " + to_type + " vcf does not contain a file")
        return False

    #Define a dict to hold the ouput
    rsid_to_pos = {}
    #Cycle through the file
    for line in islice(dbsnp_vcf_file, 1, None):
        #Convert each line to a string, remove junk characters, and split into list
        line = str(line)
        line = line.replace('\\n\'', '')
        line = line.replace('b\'', '')
        #Use if statement to rule out header lines and non-primary assembly chromosomes
        #Lines all start with b', so disregarding first two positions in each
        if line[0] != '#' and line[:2] == 'NC':
            #Split the string
            line = line.split(line_splitter_tab(line, False))
            #Check for bad trailing lines
            try:
                #Get rsid
                rsid = line[2]
                chromo = "chr" + str(int(line[0][line[0].find(".") - 2: line[0].find(".")]))
                pos = line[1]
                #Add new keys to dictionary
                rsid_to_pos[rsid] = chromo + ":" + pos
            except IndexError: #These should catch zero trailing lines
                continue
    #Verify that dictionary isn't empty
    if rsid_to_pos != {}:
        #Write dictionary to file
        with open(dict_loc, "w") as outfile:
            json.dump(rsid_to_pos, outfile)
    else:
        #Dictionary is empty, so clearly the file didn't read correctly
        print("ERROR: rsid_to_" + to_type + " creation unsuccessful.  Check input reference vcf.")
        print("ERROR: All vcf files must be in vcf format with chromosomes listed using primary assembly NCBI reference sequence, ex. chr_1 as NC_000001.10")
        return False
        
    #Output update message
    print('NOTE: rsid_to_' + to_type + ' dictionary successfully assembled and exported to ' + dict_loc)
    #Return True to indicate no problems
    return True

#Create a function that will take in the dbsnp_vcf_file_loc and make a dictionary
#for pos_to_rsid 
def make_pos_to_rsid_dict(from_type, dbsnp_vcf_file_loc, dict_loc):
    
    #Open the input file but don't read in yet
    from itertools import islice
    try:
        dbsnp_vcf_file = opener(dbsnp_vcf_file_loc)
    except FileNotFoundError:
        print("ERROR: File loc for " + from_type + " vcf does not contain a file")
        return False
    
    #Define a dict to hold the ouput
    pos_to_rsid = {}
    #Cycle through the file
    for line in islice(dbsnp_vcf_file, 1, None):
        #Convert each line to a string, remove junk characters, and split into list
        line = str(line)
        line = line.replace('\\n\'', '')
        line = line.replace('b\'', '')
        #Use if statement to rule out header lines and non-primary assembly chromosomes
        #Lines all start with b', so disregarding first two positions in each
        if line[0] != '#' and line[:2] == 'NC':
            #Split the string
            line = line.split(line_splitter_tab(line, False))
            #Check for bad trailing lines
            try:
                #Get rsid
                rsid = line[2]
                chromo = "chr" + str(int(line[0][line[0].find(".") - 2: line[0].find(".")]))
                pos = line[1]
                #Add new keys to dictionary
                pos_to_rsid[chromo + ":" + pos] = rsid
            except IndexError: #These should catch zero trailing lines
                continue
    #Verify that dictionary isn't empty
    if pos_to_rsid != {}:
        #Write dictionary to file
        with open(dict_loc, "w") as outfile:
            json.dump(pos_to_rsid, outfile)
    else:
        #Dictionary is empty, so clearly the file didn't read correctly
        print("ERROR: " + from_type + "_to_rsid creation unsuccessful.  Check input reference vcf.")
        print("ERROR: All vcf files must be in vcf format with chromosomes listed using primary assembly NCBI reference sequence, ex. chr_1 as NC_000001.10")
        return False
        
    #Output update message
    print('NOTE: ' + from_type + '_to_rsid dictionary successfully assembled and exported to ' + dict_loc)
    #Return True to indicate no problems
    return True

#Create a function that will take in two dbsnp vcf files and make a position to position converion dictionary
def make_pos_to_pos_dict(from_type, to_type, dbsnp_vcf_from_loc, dbsnp_vcf_to_loc, dict_loc):
    
    #Open the input file but don't read in yet
    from itertools import islice
    try:
        dbsnp_vcf_to = opener(dbsnp_vcf_to_loc)
    except FileNotFoundError:
        print("ERROR: File loc for " + to_type + " vcf does not contain a file")
        return False
    #Define a dict to hold the ouput
    rsid_to_pos = {}
    #Cycle through the file
    for line in islice(dbsnp_vcf_to, 1, None):
        #Convert each line to a string, remove junk characters, and split into list
        line = str(line)
        line = line.replace('\\n\'', '')
        line = line.replace('b\'', '')
        #Use if statement to rule out header lines and non-primary assembly chromosomes
        #Lines all start with b', so disregarding first two positions in each
        if line[0] != '#' and line[:2] == 'NC':
            #Split the string
            line = line.split(line_splitter_tab(line, False))
            #Check for bad trailing lines
            try:
                #Get rsid
                rsid = line[2]
                chromo = "chr" + str(int(line[0][line[0].find(".") - 2: line[0].find(".")]))
                pos = line[1]
                #Add new keys to dictionary
                rsid_to_pos[rsid] = chromo + ":" + pos
            except IndexError: #These should catch zero trailing lines
                continue
    #Verify that dictionary isn't empty
    if rsid_to_pos != {}:
        #Use rsid_to_pos dictionary for the "to" alignment to convert "from" locations to "to" locations by rsid
        #First clean-up to reduce memory demand
        del dbsnp_vcf_to
        del line
        del rsid
        del chromo
        del pos
        #Open the input file but don't read in yet
        from itertools import islice
        try:
            dbsnp_vcf_from = opener(dbsnp_vcf_from_loc)
        except FileNotFoundError:
            print("ERROR: File loc for " + from_type + " vcf does not contain a file")
            return False
    
        #Define a dict to hold the ouput
        pos_to_pos = {}
        #Cycle through the file
        for line in islice(dbsnp_vcf_from, 1, None):
            #Convert each line to a string, remove junk characters, and split into list
            line = str(line)
            line = line.replace('\\n\'', '')
            line = line.replace('b\'', '')
            #Use if statement to rule out header lines and non-primary assembly chromosomes
            #Lines all start with b', so disregarding first two positions in each
            if line[0] != '#' and line[:2] == 'NC':
                #Split the string
                line = line.split(line_splitter_tab(line, False))
                #Check for bad trailing lines
                try:
                    #Get rsid
                    rsid = line[2]
                    #Try to get "to" alignment position
                    try:
                        to_chrom_pos = rsid_to_pos[rsid]
                    except KeyError:
                        continue #Key not found, so rsid not in the to dictionary and pos_to_pos conversion not possible
                    #Get original chromosome and position
                    chromo_from = "chr" + str(int(line[0][line[0].find(".") - 2: line[0].find(".")]))
                    pos_from = line[1]
                    #Add new keys to dictionary
                    pos_to_pos[chromo_from + ":" + pos_from] = to_chrom_pos
                except IndexError: #These should catch zero trailing lines
                    continue
        #Verify that pos_to_pos dictionary isn't empty and output it if so
        if pos_to_pos != {}:
            #Write dictionary to file
            with open(dict_loc, "w") as outfile:
                json.dump(pos_to_pos, outfile)
        #Else an error has occured while reading "from" vcf
        else:
            #Dictionary is empty, so clearly the file didn't read correctly
            print("ERROR: Check" + from_type + " reference vcf. Unable to process while creating " + from_type + "_to_" + to_type + " json dictionary")
            return False
    else:
        #Dictionary is empty, so clearly the file didn't read correctly
        print("ERROR: Check" + to_type + " reference vcf. Unable to process while creating " + from_type + "_to_" + to_type + " json dictionary")
        return False
        
    #Output update message
    print('NOTE: ' + from_type + '_to_' + to_type + ' dictionary successfully assembled and exported to ' + dict_loc)
    #Return True to indicate no problems
    return True

#Make a function to open files
def opener(filename):
    f = open(filename, 'rb')
    if (f.read(2) == b'\x1f\x8b'):
        return gzip.open(filename, mode='rb')
    else:
        f.seek(0)
        return f

#Define a function that figures out how to split a line and returns the separator
def line_splitter(line, inp_separator):
    #Test different splits
    none_count = len(line.split()) 
    comma_count = len(line.split(","))
    line_count = len(line.split("|"))
    extra_slash_tab_count = len(line.split("\\t"))
    #Check input_separator if it's not false
    if inp_separator != False:
        custom_count = len(line.split(inp_separator))
    else:
        custom_count = 0 
    #Get max count value, then it's index in a dict, and then return the dict key 
    #i.e. the separator
    count_dict = {None:none_count, ",":comma_count, "|":line_count, "\\t":extra_slash_tab_count, inp_separator: custom_count}
    max_count_index = list(count_dict.values()).index(max(count_dict.values()))
    output = list(count_dict.keys())[max_count_index]
    return output      

#Define a second function that figures out how to split a line and returns the separator
#This one only considers blank spaces, tabs, and tab variants
def line_splitter_tab(line, inp_separator):
    #Test different splits
    none_count = len(line.split()) 
    extra_slash_tab_count = len(line.split("\\t"))
    #Check input_separator if it's not false
    if inp_separator != False:
        custom_count = len(line.split(inp_separator))
    else:
        custom_count = 0 
    #Get max count value, then it's index in a dict, and then return the dict key 
    #i.e. the separator
    count_dict = {None:none_count, "\\t":extra_slash_tab_count, inp_separator: custom_count}
    max_count_index = list(count_dict.values()).index(max(count_dict.values()))
    output = list(count_dict.keys())[max_count_index]
    return output 

#Define a function that gets column position and/or final sample sizes and returns a dictionary
def get_orig_col_and_sizes(trait, worked_trait_dict, align_rsid_pos, output_dir): 
    #Check if no header flag is present (assume header if no flag given)
    try: 
        no_header_flag = worked_trait_dict['no_header']
    except KeyError:
        no_header_flag = False
    #Verify/correct flag to boolean
    if type(no_header_flag) != bool:
        if type(no_header_flag) == str:
            if no_header_flag.lower() == 'true':
                no_header_flag = True
            elif no_header_flag.lower() == 'false':
                no_header_flag = False
            else:
                print("ERROR: Unrecognized no_header flag given for " + trait + ". Please select True/False.")
                return False
    #Similarly check if inp_separator was given
    try: 
        inp_separator = worked_trait_dict['separator']
    except KeyError:
        inp_separator = False
    
    #Acquire align, rsid, and pos fields from separate dict
    trait_align, rsid_name, pos_name = align_rsid_pos[trait]
    #Next try to pull non-optional field names from worked_trait_dict
    try: 
        trait_file_loc = worked_trait_dict['file_loc']
        trait_type = worked_trait_dict['type']
        chromo_name = worked_trait_dict['chromo']
        a1_name = worked_trait_dict['allele_1']
        a2_name = worked_trait_dict['allele_2']
        effect_size_name = worked_trait_dict['effect_size']
        se_name = worked_trait_dict['se']
        pval_name = worked_trait_dict['pval']
    except KeyError: 
        print("ERROR: " + trait + " config file description missing one or more of the following fields: ")
        print("ERROR: file_loc, chromo, a1, a2, effect_size, se, pval, type")
        return False
    #Read in trait_file_loc
    trait_file = opener(trait_file_loc)
    trait_raw = trait_file.readline()
    trait_file.close()
    #Call function to determine separator to use in future code
    separator = line_splitter(str(trait_raw), inp_separator)
    #Create header line to use with processing column names to numbers
    header_line = str(trait_raw)
    header_line = header_line.replace('\\n\'', '')
    header_line = header_line.replace('b\'', '')
    header_line = header_line.split(separator)
    #Create empty dictionary to take all output fields
    return_dict = {}
    #Add all fields to dict that are direct transfers from worked_trait_dict
    return_dict['inp_file'] = trait_file_loc
    return_dict['trait_type'] = trait_type
    return_dict['trait_align'] = trait_align
    return_dict['separator'] = separator
    #Process all basic fields
    return_dict['rsid_col'] = process_basic_column(rsid_name, 'rsid', no_header_flag, header_line, trait)
    return_dict['pos_col'] = process_basic_column(pos_name, 'pos', no_header_flag, header_line, trait)
    return_dict['chromo_col'] = process_basic_column(chromo_name, 'chromo', no_header_flag, header_line, trait)
    return_dict['a1_col'] = process_basic_column(a1_name, 'a1', no_header_flag, header_line, trait)
    return_dict['a2_col'] = process_basic_column(a2_name, 'a2', no_header_flag, header_line, trait)
    return_dict['effect_size_col'] = process_basic_column(effect_size_name, 'effect_size', no_header_flag, header_line, trait)
    return_dict['se_col'] = process_basic_column(se_name, 'se', no_header_flag, header_line, trait)
    return_dict['pval_col'] = process_basic_column(pval_name, 'pval', no_header_flag, header_line, trait)
    #Check if any returned -1 and return False if so (Using -1 instead of False to avoid confusion with 0s)
    if return_dict['rsid_col'] == -1 or return_dict['pos_col'] == -1 or return_dict['chromo_col'] == -1 or return_dict['a1_col'] == -1 or return_dict['a2_col'] == -1 or return_dict['effect_size_col'] == -1 or return_dict['se_col'] == -1 or return_dict['pval_col'] == -1:
        return False
    #Process the af_1 or af_2 normally too after seeing which exists
    #Read in af_1 or af_2
    if 'af_1' in worked_trait_dict.keys():
        af_1_name = worked_trait_dict['af_1']
        return_dict['af_1_col'] = process_basic_column(af_1_name, 'af_1', no_header_flag, header_line, trait)
        if return_dict['af_1_col'] < 0:
            return False
    elif 'af_2' in worked_trait_dict.keys():
        af_2_name = worked_trait_dict['af_2']
        return_dict['af_2_col'] = process_basic_column(af_2_name, 'af_2', no_header_flag, header_line, trait)
        if return_dict['af_2_col'] < 0:
            return False
    else: #Neither allele frequency given
        print("ERROR: " + trait + " config file description missing both af_1 and af_2. One must be given")
        return False
    #Process the case/control/sample_size fields
    if trait_type == 'quant':
        try:
            sample_size_name = worked_trait_dict['sample_size']
        except KeyError:
            print("ERROR: " + trait + " described as quant but sample_size lacking in config file")
            return False
    elif trait_type == 'cc':
        try:
            case_name = worked_trait_dict['cases']
            control_name = worked_trait_dict['controls']
        except KeyError:
            print("ERROR: " + trait + " described as cc but case or control counts lacking in config file")
            return False
    else: 
        print("ERROR: Unrecognized trait type selected for trait, " + trait)
        return False
    #Use header line and/or no_header_flag to convert name variables to columns
    if no_header_flag == True:
        #Depends on trait_type
        if trait_type == 'quant':
            if type(sample_size_name) != int:
                print("ERROR: No header flag set for quant trait " + trait + " but sample size is non-integer")
                return False
            #Check if sample_size_name is larger than length of header line indicating it's the number of samples
            elif sample_size_name >= len(header_line): 
                return_dict['sample_size'] = sample_size_name
            else: 
                return_dict['sample_size_col'] = sample_size_name
        else:
            #Case type must be cc at this point or it would have been caught above
            #Check if cases and controls are longer than length of header line indicating it's the number of samples
            if type(case_name) != int or type(control_name) != int:
                print("ERROR: No header flag set for cc trait " + trait + " but cases and/or controls is non-integer")
                return False
            if case_name >= len(header_line):
                return_dict['case_size'] = case_name
            else: 
                return_dict['case_col'] = case_name
            if control_name >= len(header_line):
                return_dict['control_size'] = control_name
            else:
                return_dict['control_col'] = control_name
    #Header is present in the following case, so need to map names to position in header_line
    else:
        #Depends on trait_type
        if trait_type == 'quant':
            if type(sample_size_name) == int:
                return_dict['sample_size'] = sample_size_name
            else:
                try:
                    return_dict['sample_size_col'] = header_line.index(sample_size_name)
                except ValueError:
                    print("ERROR: Given sample size for " + trait + " is non-integer and not in the header line")
                    return False
        else:
            #Case/control trait type
            #Handle cases first
            if type(case_name) == int:
                return_dict['case_size'] = case_name
            else: 
                try: 
                    return_dict['case_col'] = header_line.index(case_name)
                except ValueError:
                    print("ERROR: Given case count for " + trait + " is non-integer and not in the header line")
                    return False
            #Then handle controls
            if type(control_name) == int:
                return_dict['control_size'] = control_name
            else: 
                try: 
                    return_dict['control_col'] = header_line.index(control_name)
                except ValueError:
                    print("ERROR: Given control count for " + trait + " is non-integer and not in the header line")
                    return False
    #Make out_file (check if directory ends in slash already) and add it to return_dict
    if output_dir[-1] == '/':
        return_dict['out_file'] = output_dir + trait + "." + trait_type + ".processed.bed"
    else: 
        return_dict['out_file'] = output_dir + "/" + trait + "." + trait_type + ".processed.bed"
    #Lastly put no_header flag into the dictionary
    if no_header_flag == False:
        return_dict['no_header'] = False
    else:
        return_dict['no_header'] = True
    #Return the dictionary
    return return_dict

#Define function to process rsids and positions an get column numbers
def process_basic_column(var_name, var_type, no_header_flag, header_line, trait):
    #Denote list of key words to denote missing info
    miss_denoters = ['NA', 'None', "Null", "na", "N/a", "Na", "n/a", "Missing", "Miss", "MISS", "missing", "miss", "none"]
    #Check if no header flag is True
    if no_header_flag == True:
        #Verify that rsid and pos are both short integers and give cols if so 
        if type(var_name) == int:
            if var_name < len(header_line):
                var_col = var_name
            else: 
                print("ERROR: Given " + var_type + " column number for " + trait + " is greater than number of columns")
                return -1
        elif (var_name == None or var_name in miss_denoters) and (var_type in ['rsid', 'pos', 'chromo']):
            var_col = None
        #Appears that a non integer has been given so throw error
        else:
            print("ERROR: No-header flag selected for " + trait + ", and " + var_type + " is not a column number")
            return -1
    else:
        #Check for None
        if var_name == None or var_name in miss_denoters:
            var_col = None
        else:
            #Test the index in header_line
            try:
                var_col = header_line.index(var_name)
            except ValueError:
                print("ERROR: " + var_type + " column not in header line for " + trait)
                return -1
    return var_col

#Function that outputs a bed file for an individual trait
def out_bed_trait(trait, final_alignment, min_pval, retain_align, all_aligns, col_out_size_dict, conversion_dict):
    
    #Define a variable of the dict keys which will be used
    input_dict_keys = col_out_size_dict.keys()
    
    #############################################################################################################
    #############################################################################################################
    #############################################################################################################
    #############################################################################################################
    #We're going to create a set of if statements defining functions
    #First define all functions independent of conditions
    def a1_out(line, col_out_size_dict):
        return line[col_out_size_dict['a1_col']].upper()
    def a2_out(line, col_out_size_dict):
        return line[col_out_size_dict['a2_col']].upper()
    def effect_size_out(line, col_out_size_dict):
        return line[col_out_size_dict['effect_size_col']]
    def se_out(line, col_out_size_dict):
        return line[col_out_size_dict['se_col']]
    def pval_out(line, col_out_size_dict, min_pval):
        return max(float(line[col_out_size_dict['pval_col']]), min_pval)
    def chrom_end_out(chrom_start, a1):
        return (int(chrom_start) + len(a1))
    #Define empty variables that will be needed to call functions but may be empty
    pos_rsid_dict = None
    pos_to_pos_dict = None
    extra_aligns = None
    extra_conversion_dict = None
    self_conversion = None    
    orig_align = None
    #Define junk more_aligns_out function that otherwise won't be created if no extra alignments needed in output file
    def more_aligns_out(extra_aligns, orig_align, extra_conversion_dict, col_out_size_dict, line, self_conversion):
        return [None]
    ##############################################################################################
    #Next start conditions with pos/rsid interconversions and define chrom_start_out, chrom_out, and name_out functions
    #Use these functions to make a first_six_col_out function
    #Load dictionaries where applicable
    if col_out_size_dict['rsid_col'] == None:
        #Define functions and execute processes independent of pos to pos conversion
        #Load dictionary
        pos_to_rsid_conversion = col_out_size_dict['trait_align'] + '_to_rsid'
        with open(conversion_dict[pos_to_rsid_conversion]) as json_file:
            pos_rsid_dict = json.load(json_file)
        #Define name_out function that uses json to convert pos to rsid
        def name_out(chrom, chrom_start, pos_rsid_dict):
            return pos_rsid_dict[chrom + ":" + str(chrom_start)]
        #Check for a pos_to_pos conversion for the final alignment starting with the none case
        if col_out_size_dict['trait_align'] != final_alignment:
            #Load needed pos_to_pos dictionary
            pos_to_pos_conversion = col_out_size_dict['trait_align'] + '_to_' + final_alignment
            with open(conversion_dict[pos_to_pos_conversion]) as json_file:
                pos_to_pos_dict = json.load(json_file)
            #Define first_six_col_out function
            def first_six_col_out(line, col_out_size_dict, pos_rsid_dict = None, pos_to_pos_dict = None):
                #Get orig_chrom and orig_chrom_start
                orig_chrom = 'chr' + str(line[col_out_size_dict['chromo_col']])
                orig_chrom_start = int(line[col_out_size_dict['pos_col']])
                #Call chrom_and_chrom_start_out to get chrom and chrom_start
                chrom, chrom_start = chrom_and_chrom_start_out(orig_chrom, orig_chrom_start, pos_to_pos_dict)
                #Call name_out to get rsid
                name = name_out(orig_chrom, orig_chrom_start, pos_rsid_dict)
                #Get a1 and a2 next by calling their functions
                a1 = a1_out(line, col_out_size_dict)
                a2 = a2_out(line, col_out_size_dict)
                #Get chrom end
                chrom_end = chrom_end_out(chrom_start, a1)
                #Return list of first six columns in order
                return [chrom, chrom_start, chrom_end, name, a1, a2]
            #Define chrom_and_chrom_start_out using pos_to_pos_dict
            def chrom_and_chrom_start_out(orig_chrom, orig_chrom_start, pos_to_pos_dict):
                temp = pos_to_pos_dict[orig_chrom + ":" + str(orig_chrom_start)]
                return temp.split(":")
        #No pos to pos conversion needed to get to final_alignment
        else:
             #Define first_six_col_out function
            def first_six_col_out(line, col_out_size_dict, pos_rsid_dict = None, pos_to_pos_dict = None):
                #Call chrom_out and chrom_start_out to get chrom and chrom_start
                chrom = chrom_out(line, col_out_size_dict)
                chrom_start = chrom_start_out(line, col_out_size_dict)
                #Call name_out to get rsid
                name = name_out(chrom, chrom_start, pos_rsid_dict)
                #Get a1 and a2 next by calling their functions
                a1 = a1_out(line, col_out_size_dict)
                a2 = a2_out(line, col_out_size_dict)
                #Get chrom end
                chrom_end = chrom_end_out(chrom_start, a1)
                #Return list of first six columns in order
                return [chrom, chrom_start, chrom_end, name, a1, a2]
            #Define chrom_out and chrom_start_out simply
            def chrom_out(line, col_out_size_dict):
                return 'chr' + str(line[col_out_size_dict['chromo_col']])
            def chrom_start_out(line, col_out_size_dict):
                return int(line[col_out_size_dict['pos_col']])
    #Next is the case of rsid_to_pos interconversion
    elif col_out_size_dict['pos_col'] == None:
        #Figure out conversion dictionary for rsid_to_pos and open the dictionary
        rsid_to_pos_conversion = 'rsid_to_' + final_alignment
        with open(conversion_dict[rsid_to_pos_conversion]) as json_file:
            pos_rsid_dict = json.load(json_file)
        #Define first_six_col_out function
        def first_six_col_out(line, col_out_size_dict, rsid_to_pos_dict = None, pos_to_pos_dict = None):
            #Call name_out to get rsid
            name = name_out(line, col_out_size_dict)
            #Call chrom_and_chrom_start_out to get chrom and chrom_start from rsid
            chrom, chrom_start = chrom_and_chrom_start_out(name, rsid_to_pos_dict)
            #Get a1 and a2 next by calling their functions
            a1 = a1_out(line, col_out_size_dict)
            a2 = a2_out(line, col_out_size_dict)
            #Get chrom end
            chrom_end = chrom_end_out(chrom_start, a1)
            #Return list of first six columns in order
            return [chrom, chrom_start, chrom_end, name, a1, a2]
        #Define name_out function simply
        def name_out(line, col_out_size_dict):
            return line[col_out_size_dict['rsid_col']]
        #Define chrom_and_chrom_start_out simply
        def chrom_and_chrom_start_out(name, rsid_to_pos_dict):
            temp = rsid_to_pos_dict[name]
            return temp.split(":")
    #No pos/rsid interconversion needed
    else:
        #Check for a pos_to_pos conversion for the final alignment starting with the none case
        if col_out_size_dict['trait_align'] != final_alignment:
            #Load needed pos_to_pos dictionary
            pos_to_pos_conversion = col_out_size_dict['trait_align'] + '_to_' + final_alignment
            with open(conversion_dict[pos_to_pos_conversion]) as json_file:
                pos_to_pos_dict = json.load(json_file)
            #Define first_six_col_out function
            def first_six_col_out(line, col_out_size_dict, pos_rsid_dict = None, pos_to_pos_dict = None):
                #Get orig_chrom and orig_chrom_start
                orig_chrom = 'chr' + line[col_out_size_dict['chromo_col']]
                orig_chrom_start = int(line[col_out_size_dict['pos_col']])
                #Call chrom_and_chrom_start_out to get chrom and chrom_start
                chrom, chrom_start = chrom_and_chrom_start_out(orig_chrom, orig_chrom_start, pos_to_pos_dict)
                #Call name_out to get rsid
                name = name_out(line, col_out_size_dict)
                #Get a1 and a2 next by calling their functions
                a1 = a1_out(line, col_out_size_dict)
                a2 = a2_out(line, col_out_size_dict)
                #Get chrom end
                chrom_end = chrom_end_out(chrom_start, a1)
                #Return list of first six columns in order
                return [chrom, chrom_start, chrom_end, name, a1, a2]
            #Define name_out function simply
            def name_out(line, col_out_size_dict):
                return line[col_out_size_dict['rsid_col']]
            #Define chrom_and_chrom_start_out using pos_to_pos_dict
            def chrom_and_chrom_start_out(orig_chrom, orig_chrom_start, pos_to_pos_dict):
                temp = pos_to_pos_dict[orig_chrom + ":" + str(orig_chrom_start)]
                return temp.split(":")
        #Else
        else:
            #Define first_six_col_out function
            def first_six_col_out(line, col_out_size_dict, pos_rsid_dict = None, pos_to_pos_dict = None):
                #Call chrom_out and chrom_start_out to get chrom and chrom_start
                chrom = chrom_out(line, col_out_size_dict)
                chrom_start = chrom_start_out(line, col_out_size_dict)
                #Call name_out to get rsid
                name = name_out(line, col_out_size_dict)
                #Get a1 and a2 next by calling their functions
                a1 = a1_out(line, col_out_size_dict)
                a2 = a2_out(line, col_out_size_dict)
                #Get chrom end
                chrom_end = chrom_end_out(chrom_start, a1)
                #Return list of first six columns in order
                return [chrom, chrom_start, chrom_end, name, a1, a2]
            #Define chrom_out and chrom_start_out simply
            def chrom_out(line, col_out_size_dict):
                return 'chr' + line[col_out_size_dict['chromo_col']]
            def chrom_start_out(line, col_out_size_dict):
                return int(line[col_out_size_dict['pos_col']])
            #Define name_out function simply
            def name_out(line, col_out_size_dict):
                return line[col_out_size_dict['rsid_col']]
    #############################################################################################
    #Next we'll check for a1f vs a2f and define function a1f_out
    if 'af_1_col' in input_dict_keys:
        def a1f_out(line, col_out_size_dict):
            return [line[col_out_size_dict['af_1_col']]]
    #Then af_2_col must be present
    else:
        def a1f_out(line, col_out_size_dict):
            return [(1 - float(line[col_out_size_dict['af_2_col']]))]
    #############################################################################################
    #Check quant vs cc and then size vs. col and define the size_out function
    if col_out_size_dict['trait_type'] == 'quant':
        #Next check if sample_size or sample_size_col is in col_out_size_dict - start with actual count
        if 'sample_size' in input_dict_keys:
            def size_out(line, col_out_size_dict):
                return [col_out_size_dict['sample_size']]
        #Assume its the sample_size_col then
        else:
            def size_out(line, col_out_size_dict):
                return [line[col_out_size_dict['sample_size_col']]]
    #Else assumes trait is cc (would have been caught prior)
    else:
        #Next check whether the case_col/case_size given - start with case_size/control_size
        if 'case_size' in input_dict_keys and 'control_size' in input_dict_keys:
            def size_out(line, col_out_size_dict):
                return [col_out_size_dict['case_size'], col_out_size_dict['control_size']]
        #Assume case_size/control_col given then
        elif 'case_size' in input_dict_keys and 'control_col' in input_dict_keys:
            def size_out(line, col_out_size_dict):
                return [col_out_size_dict['case_size'], line[col_out_size_dict['control_col']]]
        #Assume case_col/control_size
        elif 'case_col' in input_dict_keys and 'control_size' in input_dict_keys:
            def size_out(line, col_out_size_dict):
                return [line[col_out_size_dict['case_col']], col_out_size_dict['control_size']]
        #Lastly case_col/control_col
        else:
            def size_out(line, col_out_size_dict):
                return [line[col_out_size_dict['case_col']], line[col_out_size_dict['control_col']]]
    ##############################################################################################
    #Lastly we're going to check for extra pos_to_pos conversions which only occurs if retain_align == True
    #Define more_aligns_out function which will need to be called inside a list
    if len(all_aligns) > 1 and retain_align == True:
        #Define the extra alignments we need to convert to
        extra_aligns = all_aligns[:]
        extra_aligns.remove(final_alignment)
        #Check if an original position was given
        if col_out_size_dict['pos_col'] != None:
            #Create orig_align variable
            orig_align = col_out_size_dict['trait_align']
            #Define a function that returns the original position  
            def get_orig_pos(line, col_out_size_dict):
                return 'chr' + str(line[col_out_size_dict['chromo_col']]) + ":" + str(line[col_out_size_dict['pos_col']])
            #Create dictionary to hold conversion dictionaries for all extra conversions
            extra_conversion_dict = {}
            #Load extra jsons by looping through extra_aligns and loading the extra dictionaries
            for extra_align in extra_aligns:
                if extra_align != orig_align:
                    #Load needed pos_to_pos dictionary
                    pos_to_pos_conversion = orig_align + '_to_' + extra_align
                    with open(conversion_dict[pos_to_pos_conversion]) as json_file:
                        extra_conversion_dict[pos_to_pos_conversion] = json.load(json_file)
            #Define self pos_to_pos conversion and its value in the dictionary (regardless of whether it's used)
            self_conversion = orig_align + '_to_' + orig_align
            extra_conversion_dict[self_conversion] = {"temp": 1}
            #Define the more_aligns_out function
            def more_aligns_out(extra_aligns, orig_align, extra_conversion_dict, col_out_size_dict, line, self_conversion):
                #Get original position
                orig_pos = get_orig_pos(line, col_out_size_dict)
                #Define list of positions to be output
                return_list = []
                #Loop through extra_aligns, get new genomic positions, and append to return_list
                for extra_align in extra_aligns:
                    #Overwrite self conversion in extra_conversion_dict (again regardless of use)
                    extra_conversion_dict[self_conversion] = {orig_pos:orig_pos}
                    #Return value using the extra_conversion_dict
                    return_list.append(extra_conversion_dict[orig_align + '_to_' + extra_align][orig_pos])
                    #Append an extra None, so we can cut it off to be consistent
                    return_list.append(None)
                #Return the list
                return return_list
        #Else we have the case of no original position, so we'll have to do all additional conversions from rsid
        else:
            #Create dictionary to hold conversion dictionaries for all extra conversions
            extra_conversion_dict = {}
            #Load extra jsons by looping through extra_aligns and loading the extra dictionaries
            for extra_align in extra_aligns:
                #Load needed rsid_to_pos dictionary
                rsid_to_pos_conversion = 'rsid_to_' + extra_align
                with open(conversion_dict[rsid_to_pos_conversion]) as json_file:
                    extra_conversion_dict[rsid_to_pos_conversion] = json.load(json_file)
            #Define the more_aligns_out function
            def more_aligns_out(extra_aligns, orig_align, extra_conversion_dict, col_out_size_dict, line, self_conversion):
                #Get rsid
                orig_rsid = line[col_out_size_dict['rsid_col']]
                #Define list of positions to be output
                return_list = []
                #Loop through extra_aligns, get new genomic positions, and append to return_list
                for extra_align in extra_aligns:
                    #Return value using the extra_conversion_dict
                    return_list.append(extra_conversion_dict['rsid_to_' + extra_align][orig_rsid])
                    #Append an extra None, so we can cut it off to be consistent
                    return_list.append(None)
                #Return the list
                return return_list
            
    #############################################################################################################
    #############################################################################################################
    #############################################################################################################
    #############################################################################################################           
    
    #Create header line list to file after checking trait_type and retain_align flag
    if col_out_size_dict['trait_type'] == 'quant':
        if retain_align == True and len(all_aligns) > 1:
            keep_aligns = all_aligns[:]
            keep_aligns.remove(final_alignment)
            header_line = ['chrom', 'chromStart', 'chromEnd', 'name', 'a1', 'a2', 'a1f', 'effect_size', 'se', 'pval', 'sample_size']
            #Cycle through keep_aligns and write out extra columns then type the return character
            for align in keep_aligns:
                header_line.append(align)
        else:
            header_line = ['chrom', 'chromStart', 'chromEnd', 'name', 'a1', 'a2', 'a1f', 'effect_size', 'se', 'pval', 'sample_size']
    else:
        if retain_align == True and len(all_aligns) > 1:
            keep_aligns = all_aligns[:]
            keep_aligns.remove(final_alignment)
            header_line = ['chrom', 'chromStart', 'chromEnd', 'name', 'a1', 'a2', 'a1f', 'effect_size', 'se', 'pval', 'cases', 'controls']
            #Cycle through keep_aligns and write out extra columns then type the return character
            for align in keep_aligns:
                header_line.append(align)
        else:
            header_line = ['chrom', 'chromStart', 'chromEnd', 'name', 'a1', 'a2', 'a1f', 'effect_size', 'se', 'pval', 'cases', 'controls']
    #Set warning flag, to let know of issues during subsequent looping
        warning_lines = 0
        dict_miss_lines = 0
    
    #Open the input file but don't read in yet
    from itertools import islice
    trait_inp_file = opener(col_out_size_dict['inp_file'])

    #Define line_count, warning lines, and dict_miss_lines
    warning_lines = 0
    dict_miss_lines = 0
    line_count = 0
    #Open write file to begin writing
    import csv
    with open(col_out_size_dict['out_file'], "w") as writeFile:
        writer = csv.writer(writeFile, delimiter = '\t')
        writer.writerow(header_line)
        #Loop through rows starting after the first and write lines to file
        for line in islice(trait_inp_file, 1, None):
            #Use try command to check for value errors due to missing data
            try:
                #Convert each line to a string, remove junk characters, and split into list
                line = str(line)
                line = line.replace('\\n\'', '')
                line = line.replace('b\'', '')
                line = line.split(col_out_size_dict['separator'])
                #Get values and store in line_out
                line_out = first_six_col_out(line, col_out_size_dict, pos_rsid_dict, pos_to_pos_dict)
                line_out.extend(a1f_out(line, col_out_size_dict))
                line_out.append(effect_size_out(line, col_out_size_dict))
                line_out.append(se_out(line, col_out_size_dict))
                line_out.append(pval_out(line, col_out_size_dict, min_pval))
                line_out.extend(size_out(line, col_out_size_dict))
                line_out.extend(more_aligns_out(extra_aligns, orig_align, extra_conversion_dict, col_out_size_dict, line, self_conversion))
                #Cut off the trailing None
                line_out.remove(None)
                #Append list to bed_file
                writer.writerow(line_out)
                #Add to line count
                line_count += 1
            except ValueError:
                warning_lines += 1
                line_count += 1
            except KeyError:
                dict_miss_lines += 1
                line_count += 1
        #Close write file
        writeFile.close()
    
    #Verify that warning_lines + dict_miss_lines < len(trait_inp_raw) i.e. at least some part of file was written
    if warning_lines + dict_miss_lines < line_count:    
        #Print update
        print(trait + " completed")
        #Check for warning message
        if warning_lines > 0:
            print("WARNING: The input file for " + trait + 
                  " had " + str(warning_lines) + " (" + "{:.1f}".format((100*warning_lines)/line_count) +"%)"
                  " line(s) with unaccepted values that were" + 
                  " excluded from the final file.  Please consider reviewing.")
        if dict_miss_lines > 0:
            print("WARNING: The input file for " + trait + 
                  " had " + str(dict_miss_lines) + " (" + "{:.1f}".format((100*dict_miss_lines)/line_count) +"%)"
                  " line(s) with rsids" + 
                  " that were missing in the given rsid-to-position conversion" + 
                  " dictionary. Please consider reviewing.")
        #Return True
        return True
    #Something must have gone wrong with the input files, so we'll output an error message
    else:
        print('ERROR: Check inputs for trait, ' + trait + '! The file was unreadable given configuration settings.') 
        return False

    ###############################################################################


#cProfile.run('main(sys.argv[1:])')
         
if __name__ == "__main__":
    main(sys.argv[1:])
