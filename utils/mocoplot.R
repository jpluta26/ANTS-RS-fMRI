
# pluta 9/20/2012
# 
# script to generate SPM style plots of 6 parameters for motion correction

# this is tested and it works, but should be added to another script- i will merge them soon
plotMocoParams <- function(mocoparamfile)
{
	mocoparams <- read.csv(mocoparamfile)


	pdf('MOCO_plot.pdf', width=12, height=8)

	# outer margins  c(bottom, left, top, right)
	par(oma=c(0,0,5,0))      


	# matrix for the order of plotting
	m <- matrix(c(1,2,3,4), 2, 2, byrow=TRUE) 

	layout(m, widths=c(12,2), heights=c(5,5), TRUE)


	n.tps <- dim(mocoparams)[1]

	# sequence for plotting
	tps <- seq(1:n.tps)   


	# 6 rigid parameters of motion correction that will be plotted
	xtrans <- mocoparams[,3]
	ytrans <- mocoparams[,4]
	ztrans <- mocoparams[,5]
	pitch  <- mocoparams[,6]
	roll   <- mocoparams[,7]
	yaw    <- mocoparams[,8]


	# setup axis parameters for translation
	ymin <- min(xtrans, ytrans, ztrans)
	ymax <- max(xtrans, ytrans, ztrans)

	# plot translation
	plot(tps, xtrans, type='l', col='blue', xlab="Time Point", ylab="Translation in mm", ylim=c(ymin, ymax))
	lines(tps, ytrans, type='l', col='green')
	lines(tps, ztrans, type='l', col='red')


	# setup a blank plot and put the legend there
	plot(c(0,1), c(0,1), ann=F, bty='n', type='n', xaxt='n', yaxt='n')
	l1 <- legend("center", xpd=TRUE, legend=c("X Trans", "Y Trans", "Z Trans"), col=c('blue', 'green', 'red'), lwd=2)

	# axis parameters for rotation
	ymin <- min(pitch, roll, yaw)
	ymax <- max(pitch, roll, yaw)


	# plot rotation
	plot(tps, pitch, type='l', col='blue', xlab="Time Point", ylab="Rotation in mm", ylim=c(ymin, ymax))
	lines(tps, roll, type='l', col='green')
	lines(tps, yaw, type='l', col='red')

	# setup another blank plot and put the second legend there
	plot(c(0,1), c(0,1), ann=F, bty='n', type='n', xaxt='n', yaxt='n')
	legend("center", xpd=TRUE, legend=c("Pitch", "Roll", "Yaw"), col=c('blue', 'green', 'red'), lwd=2)

	# figure caption
	mtext("Motion Correction Parameters", side=3, line=1, cex=2, col="black", outer=TRUE)  
	
	# write to pdf
	dev.off()
}