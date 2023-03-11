# This script creates a database that contains all the parameters that we want to run.
# For each new model, we need to provide a new functions specifying the grid.

using Pkg; Pkg.activate("../../");
using SQLite
using ProgressMeter
using ArgParse

function parse_commandline()
  s = ArgParseSettings()

  @add_arg_table! s begin
      "-m"
      arg_type = Int
      help = "Name of the model to generate scripts"
    end

  return parse_args(s)
end

function model1()
  SQLite.execute(db, """DROP TABLE IF EXISTS sourcesink1""")
  SQLite.execute(db, """
  CREATE TABLE sourcesink1 (
      beta REAL,
      gamma REAL,
      rho REAL,
      b REAL,
      cost REAL,
      mu REAL,
      PRIMARY KEY (beta, gamma, rho, b, cost, mu)
  )
  """)
  
  @showprogress for β=0.:0.02:0.2, γ= 0.9:0.1:1.1, ρ=0.1:0.15:0.40, b=0.12:0.05:0.22, c=0.:0.1:2.0
    μ = 1e-4
    params = (β, γ, ρ, b, c, μ)
    SQLite.execute(db, """INSERT INTO sourcesink1 VALUES (?, ?, ?, ?, ?, ?)""", params)
  end
end

function model2()
  SQLite.execute(db, """DROP TABLE IF EXISTS sourcesink2""")
  SQLite.execute(db, """
  CREATE TABLE sourcesink2 (
      beta REAL,
      xi REAL,
      alpha REAL,
      gamma REAL,
      rho REAL,
      eta REAL,
      b REAL,
      cost REAL,
      mu REAL,
      PRIMARY KEY (beta, xi, alpha, gamma, rho, eta, b, cost, mu)
  )
  """)
  
  @showprogress for β = 0.06:0.01:0.17, ρ = 0.005:0.005:0.1, η = 0.005:0.005:0.05, b = 0.2:0.2:1.0
    ξ, α, γ, c, μ = 1.0, 1.0, 1.0, 1e-4
    params = (β, ξ, α, γ, ρ, η, -b, c, μ)
    SQLite.execute(db, """INSERT INTO sourcesink2 VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)""", params)
  end
  
end

function model3()
  SQLite.execute(db, """DROP TABLE IF EXISTS sourcesink3""")
  SQLite.execute(db, """
  CREATE TABLE sourcesink3 (
      beta REAL,
      gamma REAL,
      rho REAL,
      b REAL,
      cost REAL,
      mu REAL,
      delta INT,
      alpha REAL,
      PRIMARY KEY (beta, gamma, rho, b, cost, mu, delta, alpha)
  )
  """)
  
  @showprogress for β=0.0:0.01:0.15, ρ=0.0:0.01:0.15, b=0.1:0.1:0.3, α=0.0:0.01:0.15
    γ = β
    c, μ, δ = 1., 0.1, 1.
    params = (β, γ, ρ, b, c, μ, δ, α)
    SQLite.execute(db, """INSERT INTO sourcesink3 VALUES (?, ?, ?, ?, ?, ?, ?, ?)""", params)
  end
end

function main()
  global db = SQLite.DB("source-sink.db")
  args = parse_commandline()
  if args["m"] == 1
    model1()
  elseif args["m"] == 2
    model2()
  elseif args["m"] == 3
    model3()
  end
end

main()
