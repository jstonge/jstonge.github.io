# teenyverse
The Teenyverse

## Workflow 

 - Make sure you have [quarto](https://quarto.org/docs/get-started/) and [julia](https://julialang.org/downloads/) installed on your machine.
 - Clone this repo.
 - Move to the directory of interests, e.g. `teenyverse/posts/source-sink` for the `source-sink` model.
 - For each model, we have the following structure:
 
```shell
├── data.json              # data file. Created when new models are ran. 
├── index.qmd              # Main file to render the interactive page. There should be only `ojs` cells in it.
├── models                 # All models live in their own directory
│   ├── source-sink1.jl
│   └── source-sink2.jl
├── sourcesink_sketch.jpg  # Sketch of the model
└── unions.jpg             # Image that is shown at the top level of the blog
```
 
### Adding new parameter configurations to an existing model
 
 - Have the model rendered somewhere, either directly on the [blog page](https://jstonge.github.io/teenyverse/posts/source-sink/) or knitting the `.qmd` document  (on VScode, this is `ctrl-shift k` when you are in the `index.qmd` file). 
 - Run new configurations as needed from the root dir of the blog page (e.g. `source-sink/`). For example, if you want to add `beta` values 0.27, `gamma` values 1.1., and the rest is the default,  you can do:
 
 ```shell
 # the rest of the params will be the default below
 julia models/source-sink1.jl --beta 0.27 -g 1.1 
 ```
 
 Running `julia models/source-sink1.jl --help` will give you the argument names and current default values:
 
 ```shell
 usage: sourcesink1.jl [--db DB] [-L L] [-O O] [--beta BETA] [-g G]
                      [-r R] [-b B] [-c C] [-m M] [-o O] [-h]

optional arguments:
  --db DB      Use Database to query parameters
  -L L         LIMIT of rows (type: Int64, default: 5)
  -O O         The OFFSET clause after LIMIT specifies how many rows
               to skip at the beginning of the result set. (type:
               Int64, default: 0)
  --beta BETA  Spreading rate from non-adopter to adopter beta (type:
               Float64, default: 0.07)
  -g G         Recovery rate gamma, i.e. rate at which adopters loose
               behavioral trait (type: Float64, default: 1.0)
  -r R         Global behavioral diffusion rho (allows the behaviour
               to spread between groups) (type: Float64, default: 0.1)
  -b B         Group benefits b (type: Float64, default: 0.18)
  -c C         Institutional cost c (type: Float64, default: 1.05)
  -m M         Noise u (type: Float64, default: 0.0001)
  -o O         Output file for results (default: ".")
  -h, --help   show this help message and exit
 ```
 Once this is done running, the rendered file should reflect the changes. To see the changes on the web page, we need to push the changes on Github. Then Github action will take care of update the web page.

### Modify/add new models

 - To modify an existing model, create a new file in `models` (see `source-sink/models/source-sink2.jl` for a template). Write the model in a function of the same name. This model should be added as a new tab in the `output` section in the main `index.qmd` file.
 - To add a new model create a new folder (see for example `teenyverse/posts/sci-group-life-cycle`). The folder should follow the same logic than first model.
