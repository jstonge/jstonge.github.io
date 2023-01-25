# This script creates a database that contains all the parameters that we want
# to run.

using Pkg; Pkg.activate("../../");
using SQLite
db = SQLite.DB("source-sink.db")

# ---------------------------------- model 1 --------------------------------- #

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

for β=0.07:0.05:0.22, γ= 0.9:0.1:1.1, ρ=0.1:0.15:0.40, b=0.12:0.05:0.22, c=.55:0.5:2.05
  p = [β, γ, ρ, b, c, 1e-4]
  params = (β, γ, ρ, b, c, 1e-4)
  SQLite.execute(db, """INSERT INTO sourcesink1 VALUES (?, ?, ?, ?, ?, ?)""", params)
end

# ---------------------------------- model 2 --------------------------------- #

SQLite.execute(db, """
CREATE TABLE sourcesink2 (
    beta REAL,
    alpha REAL,
    gamma REAL,
    rho REAL,
    b REAL,
    cost REAL,
    mu REAL,
    PRIMARY KEY (beta, alpha, gamma, rho, b, cost, mu)
)
""")

counter = 1
for β=0.01:0.01:0.25, b=0.0:0.2:1.0, c=0.5:0.5:2.0
  params = (β, 1.0, 1.0, 0.2, -b, c, 1e-4)
  SQLite.execute(db, """INSERT INTO sourcesink2 VALUES (?, ?, ?, ?, ?, ?, ?)""", params)
  counter += 1
end


# ---------------------------------- model 2 --------------------------------- #

SQLite.execute(db, """
CREATE TABLE sourcesink2 (
    beta REAL,
    alpha REAL,
    gamma REAL,
    rho REAL,
    b REAL,
    cost REAL,
    mu REAL,
    PRIMARY KEY (beta, alpha, gamma, rho, b, cost, mu)
)
""")

counter = 1
for β=0.01:0.01:0.25, b=0.0:0.2:1.0, c=0.5:0.5:2.0
  params = (β, 1.0, 1.0, 0.2, -b, c, 1e-4)
  SQLite.execute(db, """INSERT INTO sourcesink2 VALUES (?, ?, ?, ?, ?, ?, ?)""", params)
  counter += 1
end

# ---------------------------------- model 3 --------------------------------- #

#!TODO: change params and update loops
SQLite.execute(db, """
CREATE TABLE sourcesink3 (
    beta REAL,
    gamma REAL,
    rho REAL,
    b REAL,
    cost REAL,
    mu REAL,
    PRIMARY KEY (beta, gamma, rho, b, cost, mu)
)
""")

counter = 1
for β=0.01:0.03:0.61, b=0.0:0.2:1.0, c=0.5:0.5:4.0
  γ = β
  params = (β, γ, 0.2, b, c, 0.1)
  SQLite.execute(db, """INSERT INTO sourcesink3 VALUES (?, ?, ?, ?, ?, ?)""", params)
  counter += 1
end

# SQLite.execute(db, """
# DROP TABLE sourcesink2
# """)