.PHONY: clean run

clean:
		rm -f sourcesink$(model)_output/* || true
		rm -f sourcesink$(model)_output/vacc_script/* || true
		rm -f sourcesink$(model).parquet || true
		rm -f slurm-* || true

run:
		julia src/01_source-sinkDB.jl -m $(model)
		julia src/02_script2vacc.jl --db "source-sink.db" -m sourcesink$(model) -b 30
		for file in $$(ls sourcesink$(model)_output/vacc_script/*.sh); do sbatch $$file; done;
		sh check_nodes.sh
		sbatch --mem 16G --partition short --nodes 1 --ntasks=20 --time 02:59:59 --job-name=processing .run_processing.sh sourcesink$(model)_output
