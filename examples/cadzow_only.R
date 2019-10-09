x <- hlra_cadzow(1:20, r = 1)
plot(x$signal, type = "b")

y <- hlra_cadzow(list(1:20, 40:70), r = 1)
plot(y$signal[[1]], type = "b")
plot(y$signal[[2]], type = "b")


left_diags <- inv_ac_diags(11, c(.9))
right_diags <- inv_ac_diags(10, numeric(0))
z <- hlra_cadzow(1:20, r = 1, left_diags = left_diags, right_diags = right_diags)
plot(z$signal, type = "b")
