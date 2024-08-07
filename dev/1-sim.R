devtools::load_all()
p <- lsn_simpars_default(rho_scale=999, jee = 100)
print(p)

# 
s <- lsn_simulate(p)


plot(s$phi[,c("x","y")], pch = 1+s$phi$xii)
