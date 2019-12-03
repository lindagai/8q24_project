#argh worry about details like this later.

get.rv_tdt.results <- function(window.size,overlap,filepath.tped,filepath.map,results.dir){
	create.results.dir()
	get.input.files.in.windows()
	run.rv_tdt()
	}
	
	
		#i. Create directory for results
	#TODO: This could probably be a new function
	
	#ii. Read in files and get # of windows
	
	#iii. Split the map and tped files into windows of specified size
	# Read in files and get # of windows
	#Make this its own function?
	
	#iv. Get a table of snps and positions for each window
	#This has to be calculated for every method, so it might be better to have it earlier
	# in the workflow, when the number of windows is first specified?
	#NO, someone might only want to run this method without having to go through all of them
	
	#v. Run the actual RV-TDT on all functions
	
#Helper functions	
	create.results.dir()
	  get.file.param
	  n windows should be calculated as soon as the file is read in
	
	get.input.files.in.windows()
	run.rv_tdt()
	
#Helper functions
get.file.param <- function(){
	
}
	