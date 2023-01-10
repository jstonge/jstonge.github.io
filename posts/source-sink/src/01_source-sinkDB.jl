# This script creates a database that contains all the parameters that we want
# to run.

using Pkg; Pkg.activate("../../");
using SQLite
db = SQLite.DB("source-sink.db")

# ---------------------------------- model 1 --------------------------------- #

SQLite.execute(db, """
CREATE TABLE sourcesink1 (
  name REAL PRIMARY KEY,
  beta REAL,
  gamma REAL,
  rho REAL,
  b REAL,
  cost REAL,
  mu REAL,
  result TEXT 
)
""")

counter = 1
for β=0.07:0.05:0.22, γ= 0.9:0.1:1.1, ρ=0.1:0.15:0.40, b=0.12:0.05:0.22, c=.55:0.5:2.05
  p = [β, γ, ρ, b, c, 1e-4]
  params = ("sim$(counter)", β, γ, ρ, b, c, 1e-4, "sourcesink1_$(join(p, "_")).jld2")
  SQLite.execute(db, """INSERT INTO sourcesink1 VALUES (?, ?, ?, ?, ?, ?, ?, ?)""", params)
  counter += 1
end

# ---------------------------------- model 2 --------------------------------- #

# dbtest = SQLite.DB("source-sink-test.db")
# SQLite.execute(dbtest, """
# CREATE TABLE sourcesink2 (
#     name REAL PRIMARY KEY,
#     beta REAL,
#     alpha REAL,
#     gamma REAL,
#     rho REAL,
#     b REAL,
#     cost REAL,
#     mu REAL,
#     result TEXT 
# )
# """)

# counter = 1
# for α=0.5:0.1:0.6, ρ=0.1:0.15:0.40, b=0.07:0.05:0.57
#   params = ("sim$(counter)", 0.14, α, 1., ρ, b, 1.05, 1e-4, "sourcesink2_res$(counter).txt")
#   SQLite.execute(dbtest, """INSERT INTO sourcesink2 VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)""", params)
#   counter += 1
# end


SQLite.execute(db, """
CREATE TABLE sourcesink2 (
    name REAL PRIMARY KEY,
    beta REAL,
    alpha REAL,
    gamma REAL,
    rho REAL,
    b REAL,
    cost REAL,
    mu REAL 
)
""")

counter = 1
for β=0.01:0.01:0.25, b=0.0:0.2:1.0, c=0.5:0.5:2.0
  params = ("sim$(counter)", β, 1.0, 1.0, 0.2, -b, c, 1e-4)
  SQLite.execute(db, """INSERT INTO sourcesink2 VALUES (?, ?, ?, ?, ?, ?, ?, ?)""", params)
  counter += 1
end

SQLite.execute(db, """
DROP TABLE sourcesink2
""")