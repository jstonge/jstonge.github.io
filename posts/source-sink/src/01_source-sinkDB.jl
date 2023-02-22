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

# better to drop the current table before reruning

SQLite.execute(db, """
DROP TABLE sourcesink2
""")

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

# counter = 1
# for β = 0.02:0.02:0.4, ξ = 0.5:0.5:1.0, ρ = 0.0:0.1:0.4, b = 0.0:0.2:1.0
for β = 0.06:0.02:0.16, α = 1.0:1.0:2.0, ρ = 0.01:0.02:0.11, η = 0.005:0.015:0.08, b = 0.1:0.3:1.0
  # params = (β, ξ, 1.0, 1.0, ρ, 0.2, -b, 1.0, 1e-4)
  params = (β, 1.0, α, 1.0, ρ, η, -b, 1.0, 1e-4)
  SQLite.execute(db, """INSERT INTO sourcesink2 VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)""", params)
  # counter += 1
end

# ---------------------------------- model 3 --------------------------------- #

SQLite.execute(db, """
CREATE TABLE sourcesink3 (
    beta REAL,
    gamma REAL,
    rho REAL,
    b REAL,
    cost REAL,
    mu REAL,
    delta INT,
    PRIMARY KEY (beta, gamma, rho, b, cost, mu, delta)
)
""")

for β=0.02:0.02:0.4, b=0.20:0.1:1.0, ρ=0.02:0.02:0.4, δ=0:1:1
  params = (β, 0.2, ρ, b, 1.0, 0.2, δ)
  SQLite.execute(db, """INSERT INTO sourcesink3 VALUES (?, ?, ?, ?, ?, ?, ?)""", params)
end

SQLite.execute(db, """
DROP TABLE sourcesink2
""")