## Instruction to run and tweak models

#### On the VACC

The workflow to run parameter grid for a particular model is the following:

 - `src/01_source-sinkDB.jl`: This is the file where we add the parameter grid to a SQLite DB. This will be used to create bash scripts that will run on the VACC. For each model, we need to create a new section and run manually for now (i.e. run relevant function from a Julia console). In the futre we might move to a bash script using SQLite from the command line, for portability. 
 
<details>
  <summary>Click to see the code!</summary>

```julia
SQLite.execute(db, """
CREATE TABLE sourcesink2 (
    name REAL PRIMARY KEY,
    beta REAL,
    alpha REAL,
    gamma REAL,
    rho REAL,
    b REAL,
    cost REAL,
    mu REAL,
    result TEXT 
)
""")

counter = 1
for β=0.07:0.05:0.22, α=0.5:0.1:0.6, γ=0.9:0.1:1.1, ρ=0.1:0.15:0.40, b=0.12:0.05:0.22, c=.55:0.5:2.05
  params = ("sim$(counter)", β, α, γ, ρ, b, c, 1e-4, "sourcesink2_res$(counter).jld2")
  SQLite.execute(db, """INSERT INTO sourcesink2 VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)""", params)
  counter += 1
end
```
In this example, there is 864 combinations to try in our parameter sweep.
</details>

 - `src/02_script2vacc.jl`: This script is of the following form:

```shell
julia src/02_script2vacc.jl -m "sourcesink2" -b 40
```

where  `-m` is the name of the model, here `sourcesink2` (we want model name to map the filename that contains the model) and `-b` is the batch size, here `40`.

 - Now a folder should be created, in our running example it is called `sourcesink2_output`, which currently contain the bash scripts for the VACC job. All of this can be done either from the VACC or locally, but now we need to open a VACC session to call the scripts, i.e.

<details>
  <summary>Click to see the installation steps on the first time!</summary>

```shell
git clone https://github.com/jstonge/teenyverse.git
the first time we need to install relevant libraries
cd teenyverse
julia               # open the console
```

```julia
]                  # press ] to acceed the pkg management set up. You should see >Pkg instead of >Julia
instantiate        # this should install the necessery dependencies
```
</details>

```shell
ssh vacc-user1.uvm.edu -l <uvm userid>
cd teenyverse/posts/source-sink
for file in $(ls sourcesink2_output/vacc_script/*.sh); do sbatch $file; done;
# squeue -u <uvm userid>    # to check how the batch is going
```

Now the batch are sent! If everything goes according to plan, this should write the output in `sourcesink2_output/`.
