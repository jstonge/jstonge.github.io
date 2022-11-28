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
 - Move in model directory (`/models`), and run new configurations as needed. For example, if you want to add `beta` values [0.27, 0.3, 0.33], you can do:
 
 ```shell
 julia source-sink1.jl --b "0.27:0.03:0.33" # range(0.27, 0.33, step=0.03)
 ```
 
 Running `julia source-sink1.jl --help` will give you the argument names and current default values:
 
 ```shell
 usage: source-sink1.jl [--beta BETA] [--gamma GAMMA] [--rho RHO]
                       [-b B] [-c C] [-h]

optional arguments:
  --beta BETA    Spreading rate from non-adopter to adopter (type:
                 StepRangeLen, default: 0.07:0.05:0.22)
  --gamma GAMMA  Recovery rate, i.e. rate at which adopters loose
                 behavioral trait (type: StepRangeLen, default:
                 0.9:0.1:1.1)
  --rho RHO      Global behavioral diffusion (allows the behaviour to
                 spread between groups) (type: StepRangeLen, default:
                 0.1:0.15:0.4)
  -b B           Group benefits (type: StepRangeLen, default:
                 0.12:0.05:0.22)
  -c C           Institutional cost (type: StepRangeLen, default:
                 0.55:0.5:2.05)
  -h, --help     show this help message and exit
 .
 ```
 Once this is done running, the rendered file should reflect the changes. To see the changes on the web page, we need to push the changes on Github. Then Github action will take care of update the web page.

### Modify/add new models

 - To modify an existing model, create a new file in `models` (see `source-sink/models/source-sink2.jl` for a template). Write the model in a function of the same name. This model should be added as a new tab in the `output` section in the main `index.qmd` file.
 - To add a new model create a new folder (see for example `teenyverse/posts/sci-group-life-cycle`). The folder should follow the same logic than first model.
