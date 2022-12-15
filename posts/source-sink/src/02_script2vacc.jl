using Pkg; Pkg.activate("../../");
using SQLite, DataFrames, ArgParse

function parse_commandline()
    s = ArgParseSettings()
  
    @add_arg_table! s begin
        "-m"
        help = "Name of the model to generate scripts"
        "-b"
        arg_type = Int
        help = "Number of runs by batches"
      end
  
    return parse_args(s)
end
  
function main()
    args = parse_commandline()

    # MODEL_NAME = "sourcesink2"
    # NBATCHES = 80
    NBATCHES = args["b"]
    global MODEL_NAME = args["m"]
    global OUTPUT_DIR = "$(MODEL_NAME)_output"
    script_folder = "$(OUTPUT_DIR)/vacc_script"
    
    if isdir("$(OUTPUT_DIR)") == false
        mkdir("$(OUTPUT_DIR)"); mkdir("$(script_folder)") 
    end
    
    db = SQLite.DB("source-sink.db")
    c = DBInterface.execute(db, """SELECT * from $(MODEL_NAME)""") |> DataFrame

    global mem = "8gb"
    global wall_time = "02:59:59"
    global queue = parse(Int,wall_time[:2]) < 3 ? "short" : "bluemoon"
    global batch_counter = 1
    global OFFSET = 0

    function write2db(;LIMIT=0)
        open(full_script_path, "w") do io
            write(io, "#!/bin/bash\n")
            write(io, "#SBATCH --partition=$(queue)\n")
            write(io, "#SBATCH --nodes=1\n")
            write(io, "#SBATCH --mem=$(mem)\n")
            write(io, "#SBATCH --time=$(wall_time)\n")
            write(io, "#SBATCH --job-name=$(batch_counter)\n")
            write(io, "#SBATCH --output=$(OUTPUT_DIR)/res_$(batch_counter).out \n")
            write(io, "julia models/$(MODEL_NAME).jl --db 'source-sink.db' -O $(OFFSET) -L $(LIMIT) -o $(OUTPUT_DIR)")
        end
    end

    OFFSET = 0
    for i=1:nrow(c)       
        global full_script_path = "$(script_folder)/combine_folder_$(batch_counter).sh"
        if (i % NBATCHES) == 0   
            write2db(LIMIT=i)
            batch_counter += 1
            OFFSET = i
        end
    end
    write2db(LIMIT=nrow(c))
end
  
main()