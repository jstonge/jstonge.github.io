using Pkg; Pkg.activate("../../");
using SQLite, DataFrames, ArgParse

function parse_commandline()
    s = ArgParseSettings()
  
    @add_arg_table! s begin
        "--model_name"
        help = "Name of the model to generate scripts"
      end
  
    return parse_args(s)
  end
  

function main()
    args = parse_commandline()

    # MODEL_NAME = "sourcesink2"
    MODEL_NAME = args["model_name"]
    OUTPUT_DIR = "$(MODEL_NAME)_output"
    script_folder = "$(OUTPUT_DIR)/vacc_script"
    
    if isdir("$(OUTPUT_DIR)") == false
        mkdir("$(OUTPUT_DIR)"); mkdir("$(script_folder)") 
    end
    
    db = SQLite.DB("source-sink.db")
    c = DBInterface.execute(db, """SELECT * from $(MODEL_NAME)""") |> DataFrame

    mem = "8gb"
    wall_time = "02:59:59"
    queue = parse(Int,wall_time[:2]) < 3 ? "short" : "bluemoon"

    for i=1:nrow(c)
        full_script_path = "$(script_folder)/combine_folder_$(i).sh"
        open(full_script_path, "w") do io
            write(io, "#!/bin/bash\n")
            write(io, "#SBATCH --partition=$(queue)\n")
            write(io, "#SBATCH --nodes=1\n")
            write(io, "#SBATCH --mem=$(mem)\n")
            write(io, "#SBATCH --time=$(wall_time)\n")
            write(io, "#SBATCH --job-name=$(i)\n")
            write(io, "julia models/$(MODEL_NAME).jl --beta $(c[i, :beta]) -a $(c[i, :alpha]) -g $(c[i, :gamma]) -r $(c[i, :rho]) -b $(c[i, :b]) -c $(c[i, :cost]) -m $(c[i, :mu]) -o $(OUTPUT_DIR)")
        end
    end
    
end
  
main()