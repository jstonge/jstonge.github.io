.PHONY: clean run

# run until all of my nodes are done.
check_nodes_left () {
		nodes_left=`squeue -u jstonge1 | grep node | wc -l`
        while [ $nodes_left -gt 0 ]
		do
			echo 'We have $nodes_left nodes left.. sleeping'
			sleep 30
			nodes_left=`squeue -u jstonge1 | grep node | wc -l`
		done
		echo 'All nodes are done'
}

clean:
		rm $1_output/*txt
		rm $1_output/vacc_script/*
		rm $1.parquet
		rm slurm-*

run:
		julia src/01_source-sinkDB.jl -m 3
		julia src/02_script2vacc.jl --db "source-sink.db" -m "sourcesink3" -b 30
		for file in $(ls sourcesink3_output/vacc_script/*.sh); do sbatch $file; done;
		check_nodes_left
		sbatch --mem 16G --partition short --nodes 1 --ntasks=20 --time 02:59:59 --job-name=processing .run_processing.sh sourcesink3_output