eemd <-function(sig, dt, trials, nimf, noise_amp, emd_config, trials_dir=NULL)
{
	#Performs the Ensemble Empirical Mode Decomposition as described in Huang and Wu (2008) A Review on the Hilbert Huang Transform Method and its Applications to Geophysical Studies
	#It runs EMD on a given signal for N=TRIALS
 	#Each IMF set is saved to disk in TRIALS_DIR, which is created if it does not exist already.
	#Finally the EEMD function averages IMFs from all the trials together to produce an ensemble average and saves it in TRIALS_DIR 
	#INPUTS
	#	SIG is the time series to be analyzed
	#	DT is the sample rate
	#	TRIALS is the number of EMD analyses to be run
	#	NIMF is the number of IMFs to save; IMFs past this number will not be averaged
        #       NOISE_AMP determines what amplitude to make the white noise
	#	EMD_CONFIG determines how the EMD algorithm operates.
	#		EMD_CONFIG$MAX_SIFT how many times the IMFs can be sifted
	#		EMD_CONFIG$MAX_IMF maximum number of IMFs that can be generated
	#    	   	EMD_CONFIG$TOL Sifting stop criterion.
        #		EMD_CONFIG$STOP_RULE Make sure to read section on stop rules and make sure you understand what they imply!
        #		EMD_CONFIG$BOUNDARY How the start and stop of the time series are handled duing the spline process.
        #		EMD_CONFIG$SM Spline smoothing
        #		EMD_CONFIG$SPAR Smoothing parameter (only needed if sm is not none)
        #		EMD_CONFIG$WEIGHT Weight if "sm" is "spline"
	#	TRIALS_DIR is the location to store files generated during EEMD trials, if NULL then this program creates
	#	a directory called "trials" in the current directory
	#OUTPUTS are saved to TRIALS_DIR in variable EMD_RESULT
	
	if(is.null(trials_dir))
	{
		trials_dir="trials"
	}

	if(!file.exists(trials_dir))
	{
		dir.create(trials_dir)
		print(paste("Created trial directory:",trials_dir))
	}

	averaged_imfs=array(0,nimf*length(sig),dim=c(length(sig),nimf))
	averaged_noise=array(0,length(sig),dim=c(length(sig),1))
	averaged_residue=array(0,length(sig),dim=c(length(sig),1))
	for (j in 1:trials)
	{
	
		noise=runif(length(sig),min=noise_amp*-1, max=noise_amp)
		tmpsig=sig+noise
                tt=seq_len(length(sig))*dt
		emd_result=sig2imf(tmpsig,dt, emd_config)
		hht=hhtransform(emd_result)
		emd_result$hinstfreq=hht$hinstfreq
		emd_result$hamp=hht$hamp
		emd_result$noise=noise
		emd_result$original_signal=tmpsig-noise
		save(emd_result, file=paste(trials_dir, "/", "TRIAL_",sprintf("%05i",j),".RDATA",sep=""))
		if (emd_result$nimf<nimf)
		{
			trial_nimf=emd_result$nimf
		}
		else
		{
			trial_nimf=nimf
		}
		averaged_imfs[,1:trial_nimf]=averaged_imfs[,1:trial_nimf]+emd_result$imf[,1:trial_nimf]
		averaged_noise=averaged_noise+noise
		averaged_residue=averaged_residue+emd_result$residue
		print(paste("TRIAL",as.character(j),"OF",as.character(trials),"COMPLETE"))
	}
}

eemd_compile<-function(trials_dir, trials, nimf)
{
	#Averages trials together to produce a set of ensemble IMFs.
	#Produces the Hilbert spectrogram of each trial and puts it into a structure with all the other trials.
	#This is used later to generate an ensemble Hilbert spectrogram of the entire EEMD run.
	#INPUTS
	#	TRIALS_DIR is the location where the trial files produced by EEMD are stored
	#	TRIALS is the number of trials to average together
	#	NIMF is the number of IMFs to build
	#OUTPUTS
	#	EEMD_RESULT containes the ensemble IMFs


	emd_result=NULL
	if(length(list.files(trials_dir, pattern="TRIAL_\\d+\\.?RDATA$"))==0)
	{
		stop(paste("No EMD trial files found in directory", trials_dir))
	}
			
	counter=1
	sind=1
        for(file_name in list.files(trials_dir, pattern="TRIAL_\\d+\\.?RDATA$"))
        {
		load(paste(trials_dir,"/",file_name,sep=""))
		
		
                if(counter==1)
                {
			siglen=length(emd_result$original_signal)
		        averaged_imfs=array(0,nimf*siglen,dim=c(siglen,nimf))
        		averaged_noise=array(0,siglen,dim=c(siglen,1))
        		averaged_residue=array(0,siglen,dim=c(siglen,1))
			hinstfreq=array(rep(0,length(emd_result$original_signal)*nimf*trials),dim=c(length(emd_result$original_signal),nimf,trials))
			hamp=array(rep(0,length(emd_result$original_signal)*nimf*trials),dim=c(length(emd_result$original_signal),nimf,trials))
                }
                if(emd_result$nimf>=nimf)
                {
                        imf_ind=nimf
                }
                else
                {
                        imf_ind=emd_result$nimf
                }
		averaged_imfs[,1:imf_ind]=averaged_imfs[,1:imf_ind]+emd_result$imf[,1:imf_ind]
		averaged_noise=averaged_noise+emd_result$noise
		averaged_residue=averaged_residue+emd_result$residue

                hinstfreq[,1:imf_ind,counter]=emd_result$hinstfreq[,1:imf_ind]
                hamp[,1:imf_ind,counter]=emd_result$hamp[,1:imf_ind]
                sind=sind+imf_ind

                counter=counter+1
                if(counter>trials)
                {
                        break
                }
        }
	counter=counter-1
	averaged_imfs=averaged_imfs/counter
	averaged_noise=averaged_noise/counter
	averaged_residue=averaged_residue/counter
	
        if(counter<trials)
        {
                warning("Number of trials requested was greater than the number of trials found in the trials directory")
        }

	eemd_result=c()
        eemd_result$nimf=nimf
        eemd_result$dt=emd_result$dt
        eemd_result$original_signal=emd_result$original_signal
	eemd_result$averaged_imfs=averaged_imfs
	eemd_result$averaged_noise=averaged_noise
	eemd_result$averaged_residue=averaged_residue
        eemd_result$hinstfreq=hinstfreq
        eemd_result$hamp=hamp
	eemd_result$trials=trials 
        invisible(eemd_result)
}

eemd_resift <- function(eemd_result, emd_config, resift_rule)
{
	#Resifts averaged IMFs generated by EEMD to generate valid IMFs for Hilbert Transform
	#INPUTS
	#	EEMD_RESULT contains the averaged IMF set generated from EEMD.
	#	EMD_PARAMS controls how the EMD of the averaged IMFs is handled.
        #         EMD_CONFIG$MAX_SIFT how many times the IMFs can be sifted
        #         EMD_CONFIG$MAX_IMF maximum number of IMFs that can be generated
        #         EMD_CONFIG$TOL Sifting stop criterion.
        #         EMD_CONFIG$STOP_RULE Make sure to read section on stop rules and make sure you understand what they imply!
        #         EMD_CONFIG$BOUNDARY How the start and stop of the time series are handled duing the spline process.
        #         EMD_CONFIG$SM Spline smoothing
        #         EMD_CONFIG$SPAR Smoothing parameter (only needed if sm is not none)
        #         EMD_CONFIG$WEIGHT Weight if "sm" is "spline"
	#	RESIFT_RULE determines how the resifting occurs
	#	If resift_rule is numeric, get the nth IMF (so if resift_rule is 2, get the 2nd resifted IMF)
	#	If resift_rule is "last", return the last IMF.
	#	If resift_rule is "max_var", return the IMF with the most variance
	#	If resift_rule is "all", get all the IMFs returned by rerunning EMD on the averaged IMFs made by EEMD.
	#	This will likely be quite large.
	#OUTPUTS
	#	EEMD_RESULT$IMF a set of IMFs that are generated from the EEMD imfs.

	resift_result=eemd_result
	resift_result$imf=c()
	resift_result$averaged_imfs=NULL

	if(!is.numeric(resift_rule) & !resift_rule %in% c("last", "max_var", "all"))
        {
               e=simpleError(paste("Did not recognize resift_rule:", resift_rule))
               stop(e)
        }		
	for(k in seq_len(dim(eemd_result$averaged_imfs)[2]))
	{
		if(sum(eemd_result$averaged_imfs[,k]==0)!=length(eemd_result$averaged_imfs[,k]))
		{
			emd_result=sig2imf(eemd_result$averaged_imfs[,k], eemd_result$dt, emd_config)
			if(is.numeric(resift_rule))
			{
				if(emd_result$nimf>=resift_rule)
				{
					resift_result$imf=cbind(resift_result$imf, emd_result$imf[,resift_rule])
				}
				else
				{
					resift_result$imf=cbind(resift_result$imf, NA)
				}
			}
			else
			{
				if(resift_rule=="last")
				{
					resift_result$imf=cbind(resift_result$imf, emd_result$imf[,emd_result$nimf])
				}
				
				if(resift_rule=="max_var")
				{
					var_list=c()
					for(j in seq_len(emd_result$nimf))
					{
						var_list=cbind(var_list, var(emd_result$imf[,j]))
					}
					
					resift_result$imf=cbind(resift_result$imf, emd_result$imf[,var_list==max(var_list)])
				}

				if(resift_rule=="all")
				{
					resift_result$imf=cbind(resift_result$imf, emd_result$imf)
				}
			}	
		}
	}
	
	resift_result$resift_emd_config=emd_config
	resift_result$resift_rule=resift_rule
	resift_result$nimf=dim(resift_result$imf)[2]
	resift_result=hhtransform(resift_result)
	invisible(resift_result)
}
		

hh_render <- function(hres,max_freq, freq_step, imf_list=NULL)
{
	#Renders a spectrogram of EMD or Ensemble EMD (EEMD) results.
	#INPUTS
	#	HRES is a matrix of data generated by EEMD_COMPILE or the output of HHTRANSFORM
	#	it represents a set on all time/frequency/amplitude points from the given EEMD run
        #       MAX_FREQ is the maximum frequency span to plot.
	#	FREQ_STEP is the frequency resolution
        #       IMF_LIST controls which IMFs to include in the spectrogram
	#OUTPUTS
	#	HSPEC is a spectrogram matrix ready to be plotted by HHSPEC_IMAGE
	#Danny Bowman
	#UNC Chapel Hill
	#August 2012

        if(is.null(imf_list))
        {
            imf_list=1:hres$nimf
        }

	if(!"trials" %in% names(hres))
 	{
		hres$trials=1
	}

	#If the spectra is generated from EMD instead of EEMD, make it resemble EEMD but with only 1 trial.
	if(length(dim(hres$hamp))==2)
	{
		hres$hamp=array(hres$hamp, dim=c(dim(hres$hamp)[1], dim(hres$hamp)[2], 1))
		hres$hinstfreq=array(hres$hinstfreq, dim=c(dim(hres$hamp)[1], dim(hres$hamp)[2], 1))
	}

	x_len=-1
	counter=1
	hspec=c()
	x_len=dim(hres$hamp)[1]
	y_len=max_freq/freq_step
	hspec$z=array(rep(0,(x_len*y_len)),dim=c(x_len,y_len))
	hspec$cluster=array(rep(0,(x_len*y_len)),dim=c(x_len,y_len)) #Shows how many times a given grid node has data.
	r=1
	s=0
	print("Rendering spectrogram...")
	percentage=max_freq/100
	print_percent=percentage*10
        freqs=array(hres$hinstfreq[,imf_list,], dim=c(length(hres$original_signal),length(imf_list)*hres$trials)) #Get only the requested IMFs
        amps=array(hres$hamp[,imf_list,], dim=c(length(hres$original_signal),length(imf_list)*hres$trials)) #Get only the requested IMFs
	repeat
	{
		atmp=amps
		atmp[!(freqs>s & freqs<=s+freq_step)]=0 #Set all amplitudes not in frequency range to 0
		hspec$z[,r]=rowSums(atmp) #Add them to image
		trial_counter=0
		cluster=array(rep(0,hres$trials*x_len),dim=c(x_len,hres$trials))
		for(k in seq_len(hres$trials))
		{
                        if(length(imf_list)>1)
                        {
			    cluster[,k]=rowSums(atmp[,((k-1)*length(imf_list)+1):(k*length(imf_list))])
                        }
                        else
                        {
                            cluster[,k]=atmp[,((k-1)*length(imf_list)+1):(k*length(imf_list))]
                        }
			cluster[cluster[,k]>0,k]=1 #Set all amplitudes in frequency to 1
	}
		hspec$cluster[,r]=rowSums(cluster) #Add them together to give a sum total of trials that yield a signal in a given grid cell
		r=r+1
		s=s+freq_step
		if(s>=max_freq-freq_step)
		{
			break
		}	
		if(print_percent<s)
		{
			print(paste(as.character(print_percent/percentage)," percent complete..."))
			print_percent=print_percent+percentage*10
		}
	}
	
	hspec$x=seq_len(x_len)*hres$dt
	hspec$y=seq_len(y_len)*freq_step
	hspec$z=hspec$z/counter
	hspec$trials=counter
	hspec$original_signal=hres$original_signal
	hspec$max_freq=max_freq
	hspec$freq_step=freq_step
	hspec$dt=hres$dt
	invisible(hspec) #Return the spectrogram structure.
}

ftspec_image <- function(ft, time_span, freq_span, amp_span, amp_units=NULL, amp_unit_conversion=NULL, grid=TRUE, colorbar=TRUE, backcol=c(0, 0, 0), pretty=TRUE, cex=1, main="")
{
	#Plots a Fourier spectrogram
	#INPUTS
	#FT is the Fourier transform input parameters, adopted from Jonathan Lees' code in RSEIS
	#	FT$XT is the signal
	#	FT$DT is the sample rate
	#	FT$NFFT is the fft length
	#	FT$NS is the number of samples in a window
	#	FT$NOV is the number of samples to overlap
        #TIME_SPAN is the time span to plot, [0,-1] plots everything
        #FREQ_SPAN is the frequency span to plot (<=max frequency in spectrogram), [0,-1] plots everything up to the Nyquist frequency
	#AMP_SPAN is the amplitude range to plot.  [0, -1] plots everything.
        #AMP_UNITS is the units of amplitude, used for axes labels
        #AMP_UNIT_CONVERSION tells how to convert units in input data to units to display
        #GRID is a boolean asking whether to display grid lines
        #COLORBAR is a boolean asking whether to plot an amplitude colorbar
        #BACKCOL is a 3 element vector of RGB values for the background of the spectrogram, based on a 0 to 255 scale: [red, green, blue]
        #PRETTY is a boolean asking whether to adjust axis labels so that they're pretty (TRUE) or give the exactly specified time and frequency intervals (FALSE)
        #CEX is the size of the text on the figure
        #MAIN is the title of the plot.
	
	
	options(digits.secs=3)

        if(time_span[2]<0)
        {
                time_span[2]=length(ft$xt)*ft$dt
        }

        if(time_span[2]>length(ft$xt)*ft$dt)
        {
                time_span[2]=length(ft$xt)*ft$dt
                warning("Requested time window is longer than the actual signal.")
        }

        if(freq_span[2]<0)
        {
                freq_span[2]=(1/ft$dt)/2
        }
        if(freq_span[2]>(1/ft$dt)/2)
        {
                freq_span[2]=(1/ft$dt)/2
                warning("Requested frequency window is higher than the Nyquist frequency.")
        }

        if(is.null(amp_unit_conversion))
        {
                amp_unit_conversion=1
        }

        if(pretty)
        {
            time_labels=pretty(seq(0,time_span[2]-time_span[1],length.out=10)+time_span[1], n=10)
            freq_labels=pretty(seq(freq_span[1],freq_span[2],length.out=5), n=5)
            time_span=c(min(time_labels), max(time_labels))
            freq_span=c(min(freq_labels), max(freq_labels))
        }
        else
        {
           time_labels=format(seq(0,time_span[2]-time_span[1],length.out=10)+time_span[1],digits=2)
           freq_labels=format(seq(freq_span[1],freq_span[2],length.out=5), digits=2)
        }

        tind=seq_len((time_span[2]-time_span[1])/ft$dt)+(time_span[1])/ft$dt

        par(mai=c(0.1, 0.1, 0.1, 0.1))
        plot(c(-0.15,1),c(-0.15,1),type="n",axes=FALSE,xlab="", ylab="") # Set up main plot window
        par(mai=c(0.1, 0.1, 0.1, 0.1))

        xt=ft$xt[tind]*amp_unit_conversion


	ev=evolfft(xt, dt=ft$dt, Nfft=ft$nfft, Ns=ft$ns, Nov=ft$nov, fl=freq_span[1], fh=freq_span[2])

        ampcex=0.75

        #Plot TRACE

        xt=ft$xt[tind]*amp_unit_conversion
        tt=tind*ft$dt-time_span[1]
        trace_y=0.75
        trace_x=0
        trace_yspan=0.10
        trace_xspan=0.9
        trace_at=seq(trace_y,trace_y+trace_yspan,length.out=2)
        trace_labels=format(seq(min(xt), max(xt), length.out=2), digits=2)
        trace_scale=trace_yspan/(max(xt)-min(xt))
        tt_scale=trace_xspan/max(tt)
        axis(4,pos=trace_x+trace_xspan,at=trace_at, labels=trace_labels, cex.axis=ampcex)
        lines((tt*tt_scale+trace_x), (trace_y+(xt-min(xt))*trace_scale))
        rect(trace_x, trace_y, trace_x+trace_xspan, trace_y+trace_yspan)

        #Plot IMAGE
	
        f_ind=(ev$freqs>=freq_span[1] & ev$freqs <= freq_span[2])
        freqs=ev$freqs[f_ind]
        image_z=t(ev$DSPEC[1:(ev$numfreqs/2),])
        image_z=image_z[,f_ind]
        if(amp_span[2]<0)
        {

             amp_span[1]=min(image_z)
             amp_span[2]=max(image_z)
        }
        image_z[image_z<amp_span[1]]=NA
        image_z[image_z>amp_span[2]]=amp_span[2]
        image_y=0
        image_x=0
        image_yspan=0.75
        image_xspan=0.9
        image_xvec=seq(image_x, image_x+image_xspan, length.out=length(ev$tims))
        image_yvec=seq(image_y, image_y+image_yspan, length.out=length(freqs))
        time_at=seq(image_x,image_x+image_xspan,length.out=length(time_labels))
        freq_at=seq(image_y,image_y+image_yspan, length.out=length(freq_labels))
        colorbins=500
        colormap=rainbow(colorbins,start=0,end=5/6)
        rect(image_x,image_y,image_x+image_xspan,image_y+image_yspan,col=rgb(red=backcol[1], green=backcol[2], blue=backcol[3], maxColorValue=255))
        image(image_xvec,image_yvec, image_z, col=colormap,add=TRUE)
        axis(2, pos=image_x, at=freq_at,labels=freq_labels, cex.axis=cex)
        axis(1, pos=image_y, at=time_at,labels=time_labels, cex.axis=cex)
        
        #Plot GRID
        
        if(grid)
        {
                line_color=rgb(red=100, green=100, blue=100, maxColorValue=255)
                line_type=3
                for(k in 2:(length(time_at)-1))
                {
                        lines(c(time_at[k], time_at[k]), c(image_y, trace_y+trace_yspan), col=line_color, lty=line_type)
                }

                for(k in 2:(length(freq_at)-1))
                {
                        lines(c(image_x, image_x+image_xspan), c(freq_at[k], freq_at[k]), col=line_color, lty=line_type)
                }
        }

        #Plot COLORBAR

        if(colorbar)
        {
                color_x=image_x+image_xspan+0.005
                color_xspan=0.025
                color_y=image_y+image_yspan-0.20
                color_yspan=0.10
                color_xvec=c(color_x,color_x+color_xspan)
                color_yvec=seq(color_y, color_y+color_yspan, length.out=colorbins)
                color_at=seq(color_y,color_y+color_yspan,length.out=2)
                color_labels=format(seq(amp_span[1], amp_span[2],length.out=2),digits=2)
                colorbar=array(seq_len(colorbins),dim=c(1, colorbins))
                image(color_xvec, color_yvec, colorbar, col=colormap, axes=FALSE, add=TRUE)
                #axis(4,pos=color_x+color_xspan, at=color_at, labels=color_labels, cex.axis=ampcex)
        }

        #Plot TEXT
	
	text(image_x-0.095, image_y+image_yspan/2, srt=90, "Frequency", cex=cex)
        text(image_x+image_xspan/2, image_y-0.1, "Time", cex=cex)
        text(trace_x+trace_xspan/2, trace_y+trace_yspan+0.05,main, cex=cex)
        text(color_x+0.0125, color_y-0.0125, format(amp_span[1], digits=2), cex=ampcex)
        text(color_x+0.0125, color_y+color_yspan+0.0125, format(amp_span[2], digits=2), cex=ampcex)
        if(!is.null(amp_units))
        {
                text(trace_x+trace_xspan/2, trace_y+trace_yspan+0.025,paste("Trace and Spectrogram Amplitudes in",amp_units), cex=ampcex)
        }
 
       #Plot WINDOW
       
       rwidth=trace_xspan*(ev$wpars$Ns/length(xt))		
       rect(trace_x+trace_xspan-rwidth, trace_y+trace_yspan, trace_x+trace_xspan, trace_y+trace_yspan+0.01, col="black")

}
	
hhspec_image <- function(hspec,time_span,freq_span, amp_span,cluster_span=NULL, amp_units=NULL, amp_unit_conversion=NULL, grid=TRUE, colorbar=TRUE, backcol=c(0, 0, 0), pretty=TRUE, cex=1, main="")
{
	#Plots a spectrogram of the EEMD processed signal as an image.	
	#INPUTS
	#	HSPEC is the subsetted spectrogram  from HH_RENDER.
	#		HSPEC$X is time
	#		HSPEC$Y is frequency
	#		HSPEC$Z is amplitude normalized to trials
	#		HSPEC$CLUSTER is a matrix containing integer values corresponding to the number of times a signal was recorded in a given spectrogram cell during EEMD
	#		The more often the signal is recorded, the more likely it is that the signal is real and not noise
	#		HSPEC$TRIALS is the number of times EEMD was run to generate signal
	#		HSPEC$ORIGINAL_SIGNAL is the original seismogram (without added noise)
	#		HSPEC$MAX_FREQ is the top of the plotted frequency range
	#		HSPEC$FREQ_STEP is the frequency discretization of the spectrogram
	#		HSPEC$DT is the sample rate
	#	TIME_SPAN is the time span to plot, [0,-1] plots everything
	#	FREQ_SPAN is the frequency span to plot (<=max frequency in spectrogram), [0,-1] plots everything
	#	AMP_SPAN is the amplitude span to plot, everything below is set to black, everything above is set to max color, [0, -1] scales to range in signal
	#	CLUSTER_SPAN plots only the parts of the signal that have a certain number of data points per pixel [AT LEAST, AT MOST] this only applies to EEMD with multiple trials.
	#	AMP_UNITS is the units of amplitude, used for axes labels
	#	AMP_UNIT_CONVERSION tells how to convert units in input data to units to display
	#	GRID is a boolean asking whether to display grid lines
	#	COLORBAR is a boolean asking whether to plot an amplitude colorbar
        #       BACKCOL is a 3 element vector of RGB values for the background of the spectrogram, based on a 0 to 255 scale: [red, green, blue]
        #       PRETTY is a boolean asking whether to adjust axis labels so that they're pretty (TRUE) or give the exactly specified time and frequency intervals (FALSE)
        #	CEX is character size for labels and titles
        #	SUBFIGTEXT is a hack to get subfigure labels.
	#
	#OUTPUTS
	#	HHIMAGE is the result of the requested spectrogram subsetting.


        options(digits.secs=3)
	
	if(time_span[2]<0)
	{
		time_span[2]=length(hspec$original_signal)*hspec$dt
	}
	
	if(time_span[2]>length(hspec$original_signal)*hspec$dt)
	{
		time_span[2]=length(hspec$original_signal)*hspec$dt
		warning("Requested time window is longer than the actual signal.")
	}
	
	if(freq_span[2]<0)
	{
		freq_span[2]=hspec$max_freq
	}
	if(freq_span[2]>hspec$max_freq)
	{
		freq_span[2]=hspec$max_freq
		warning("Requested frequency window is higher than maximum frequency in the spectrogram.")
	}

	if(is.null(amp_unit_conversion))
	{
		amp_unit_conversion=1
	}
	
	
	if(pretty)
        {
            time_labels=pretty(seq(0,time_span[2]-time_span[1],length.out=10)+time_span[1], n=10)
            freq_labels=pretty(seq(freq_span[1],freq_span[2],length.out=5), n=5)
            time_span=c(min(time_labels), max(time_labels))
            freq_span=c(min(freq_labels), max(freq_labels))
        }
        else
        {
           time_labels=format(seq(0,time_span[2]-time_span[1],length.out=10)+time_span[1],digits=2)
           freq_labels=format(seq(freq_span[1],freq_span[2],length.out=5), digits=2)
        }

	f_ind=seq_len((freq_span[2]-freq_span[1])/hspec$freq_step)+(freq_span[1]/hspec$freq_step)
        tind=seq_len((time_span[2]-time_span[1])/hspec$dt)+(time_span[1])/hspec$dt
	
	clustermatrix=hspec$cluster[tind,f_ind]
	image_z=hspec$z[tind,f_ind]*amp_unit_conversion
	if(!is.null(cluster_span))
	{
		image_z[clustermatrix<cluster_span[1] | clustermatrix>cluster_span[2]]=NaN
	}
	image_z=image_z/clustermatrix #Scale by number of times bin got values
	if(amp_span[2]<0)
	{
		amp_span[2]=max(image_z[!is.nan(image_z)])
	}
	image_z[image_z<amp_span[1]]=NA
	image_z[image_z>amp_span[2]]=amp_span[2]
	image_z[1,1]=amp_span[2]
	maintitle=paste(hspec$network,"Station",hspec$STN,"Component",hspec$COMP,"Original Signal")
	par(mai=c(0.1, 0.1, 0.1, 0.1))
	plot(c(-0.15,1),c(-0.15,1),type="n",axes=FALSE,xlab="", ylab="") # Set up main plot window
	par(mai=c(0.1, 0.1, 0.1, 0.1))
	
	ampcex=0.75

	#Plot TRACE
	
        xt=hspec$original_signal[tind]*amp_unit_conversion
        tt=tind*hspec$dt-time_span[1]
	trace_y=0.75
	trace_x=0
	trace_yspan=0.10
	trace_xspan=0.9
	trace_at=seq(trace_y,trace_y+trace_yspan,length.out=2)
	trace_labels=format(seq(min(xt), max(xt), length.out=2), digits=2)
	trace_scale=trace_yspan/(max(xt)-min(xt))
	tt_scale=trace_xspan/max(tt)
	axis(4,pos=trace_x+trace_xspan,at=trace_at, labels=trace_labels, cex.axis=ampcex)
	lines((tt*tt_scale+trace_x), (trace_y+(xt-min(xt))*trace_scale))
	rect(trace_x, trace_y, trace_x+trace_xspan, trace_y+trace_yspan)

	#Plot IMAGE

	image_y=0
	image_x=0
	image_yspan=0.75
	image_xspan=0.9
	image_xvec=seq(image_x, image_x+image_xspan, length.out=length(tind))
	image_yvec=seq(image_y, image_y+image_yspan, length.out=length(f_ind))
	time_at=seq(image_x,image_x+image_xspan,length.out=length(time_labels))
	freq_at=seq(image_y,image_y+image_yspan, length.out=length(freq_labels))
	colorbins=500
	colormap=rainbow(colorbins,start=0,end=5/6)
	rect(image_x,image_y,image_x+image_xspan,image_y+image_yspan,col=rgb(red=backcol[1], green=backcol[2], blue=backcol[3], maxColorValue=255))
	image(image_xvec,image_yvec, image_z, col=colormap,add=TRUE)
	axis(2, pos=image_x, at=freq_at,labels=freq_labels, cex.axis=cex)
	axis(1,pos=image_y, at=time_at,labels=time_labels, cex.axis=cex)


        #Plot GRID
	if(grid)
	{
        	line_color=rgb(red=100, green=100, blue=100, maxColorValue=255)
		line_type=3
        	for(k in 2:(length(time_at)-1))
        	{
                	lines(c(time_at[k], time_at[k]), c(image_y, trace_y+trace_yspan), col=line_color, lty=line_type)
        	}
        
        	for(k in 2:(length(freq_at)-1))
        	{
                	lines(c(image_x, image_x+image_xspan), c(freq_at[k], freq_at[k]), col=line_color, lty=line_type)
        	}
	}

	#Plot COLORBAR
	
	if(colorbar)
	{
		color_x=image_x+image_xspan+0.005
        	color_xspan=0.025
		color_y=image_y+image_yspan-0.20
        	color_yspan=0.10
		color_xvec=c(color_x,color_x+color_xspan)
		color_yvec=seq(color_y, color_y+color_yspan, length.out=colorbins)
		color_at=seq(color_y,color_y+color_yspan,length.out=2)
		color_labels=format(seq(amp_span[1],amp_span[2],length.out=2),digits=2)
		colorbar=array(seq_len(colorbins),dim=c(1, colorbins))
        	image(color_xvec, color_yvec, colorbar, col=colormap, axes=FALSE, add=TRUE)
	}
	
	#Plot TEXT
	#text(color_x+color_xspan+0.05,color_y+color_yspan/2,srt=90,paste("Spectrum Amplitude (", amp_units, ")", sep=""))
	#text(trace_x+trace_xspan+0.05,trace_y+trace_yspan/2,srt=90,paste("Trace Amplitude (", amp_units, ")", sep=""))
	text(image_x-0.095, image_y+image_yspan/2, srt=90, "Frequency", cex=cex)
	text(image_x+image_xspan/2, image_y-0.1, "Time", cex=cex)
	text(trace_x+trace_xspan/2, trace_y+trace_yspan+0.05,main, cex=cex)
	text(color_x+0.0125, color_y-0.0125, format(amp_span[1], digits=2), cex=ampcex)
	text(color_x+0.0125, color_y+color_yspan+0.0125, format(amp_span[2], digits=2), cex=ampcex)
	if(!is.null(amp_units))
	{
		text(trace_x+trace_xspan/2, trace_y+trace_yspan+0.025,paste("Trace and Spectrogram Amplitudes in",amp_units), cex=ampcex)
	}
	
}
	
hhtransform <- function(imf_set)

{
    #Transform IMFs into instantaneous frequency and amplitude using the Hilbert transform
    #INPUTS
    #EMD_RESULT is the EMD decomposition of a signal invisibleed by sig2imf.R
    #Danny Bowman
    #OUTPUTS
    #  HHT is the Hilbert Transform of EMD_RESULT
    #UNC Chapel Hill
   
    
   if("averaged_imfs" %in% names(imf_set))
    {
        imf_set$imf=imf_set$averaged_imfs
    }

    tt=seq(from=0, by=imf_set$dt, length=length(imf_set$original_signal))
    hht_result=imf_set
    hilbert=hilbertspec(imf_set$imf,tt=tt)
    hht_result$hamp=hilbert$amplitude
    hht_result$hinstfreq=hilbert$instantfreq
    invisible(hht_result)
}

plot_imfs <-function(sig, time_span, imf_list, original_signal, residue, fit_line=FALSE, lwd=1, cex=1)
{
    #Better IMF plotter
    #This function plots IMFs on the same plot so they can be checked for mode mixing or other problems.
    #It plots shifted traces in a single window
    #INPUTS
    #    SIG is the signal data structure returned by EEMD or EMD analysis
    #    Note that SIG$AVERAGED_IMFS will be plotted instead of SIG$IMF, likewise SIG$AVERAGED_RESIDUE takes precidence
    #    over SIG$RESIDUE, if both exist.
    #        SIG$IMF is a N by M array where N is the length of the signal and M is the number of IMFs
    #        SIG$ORIGINAL_SIGNAL is the original signal before EEMD
    #        SIG$RESIDUE is the residual after EMD
    #        SIG$DT is the sample rate
    #    TIME_SPAN is a 2 element vector giving the time range to plot
    #    IMF_LIST is the IMFs to plot
    #    ORIGINAL_SIGNAL is a boolean asking if you are going to plot the original signal also (defaults to be on top)
    #    RESIDUE is a boolean asking if you are going to plot the residual (defaults to be on bottom)
    #	 FIT_LINE is a boolean asking if you want to plot a line showing the sum of IMFs on top of the original signal (to check how the selected IMFs reconstruct the original signal)
    #	 LWT is the line weight (for plotting figures)
    #    CEX is the size of text (for plotting figures)
   
 
    if(time_span[2]<0)
    {
        time_span[2]=length(sig$original_signal)*sig$dt
    }
    
    if(time_span[1]==0)
    {
        time_span[1]=sig$dt
    }
    
    if("averaged_imfs" %in% names(sig))
    {
        sig$imf=sig$averaged_imfs
    }

    if("averaged_residue" %in% names(sig))
    {
        sig$residue=sig$averaged_residue
    }


    time_ind=ceiling(time_span[1]/sig$dt):floor(time_span[2]/sig$dt)
    tt=time_ind*sig$dt
    
    plot(c(0, 1), c(0, 1), type="n", axes=FALSE, xlab="Time (s)", ylab="", cex.lab=cex)
    
    yax_labs=c()
    snum=length(imf_list)+residue+original_signal
    sp=1/snum # Spacing of subplots

    if(original_signal)
    {
         snum=snum+1
         os=sig$original_signal[time_ind]-mean(sig$original_signal[time_ind])
         scale=max(abs(os)) 
    }
    else
    {
        scale=max(abs(sig$imf))
    }
    
    if(residue)
    {
        snum=snum+1
        res=sig$residue[time_ind]-mean(sig$residue[time_ind])
	res=res*(sp/(2*scale))
        yax_labs=append(yax_labs,"Residue")
    }

    
    trace_pos=sp/2 #Where the trace starts on the plot
    imfs=sig$imf*(sp/(scale*2))
    ts=(tt-min(tt))*(1/(time_span[2]-time_span[1]))

    if(residue)
    {
        lines(ts, res+trace_pos, lwd=lwd)
        trace_pos=trace_pos+sp
    } 
    
    for(k in rev(imf_list))
    {
       lines(ts, imfs[time_ind,k]+trace_pos, lwd=lwd)
       trace_pos=trace_pos+sp
       yax_labs=append(yax_labs, paste("IMF",k))
    }
    if(original_signal)
    {
        lines(ts, os*(sp/(2*scale))+trace_pos, lwd=lwd)
        yax_labs=append(yax_labs,"Signal")
        if(fit_line)
        {
            if(length(imf_list)>1)
            {
                fline=rowSums(imfs[time_ind,imf_list])
            }
            else
            {
                fline=imfs[time_ind,imf_list]
            }
            if(residue)
            {
                fline=fline+res
            }
            lines(ts, fline+trace_pos, lwd=lwd, col="red")
        }
    }
    xax_labs=pretty(seq(min(tt)-sig$dt, max(tt), length.out=11))
    axis(1, pos=0, at=seq(0,1, length.out=length(xax_labs)), labels=xax_labs, cex.axis=cex)
    axis(2, pos=0, at=seq(sp/2, 1, by=sp), labels=yax_labs, cex.axis=cex)
    segments(c(0,0,1, 0), c(0, 1, 1, 0), c(0,1, 1, 1), c(1,1, 0, 0), lwd=lwd) 
}

sig2imf <- function(sig, dt, emd_config)

{
    #Extract IMFs
    #This function is intended to take data, recover IMFs, and save them
    #It calls and runs code developed by Donghoh Kim, Hee-Seok Oh as part of the "EMD" package available on CRAN.
    #I have modified their code and included some of it in this repository.
    #INPUTS
    #	SIG is the time series
    #   DT is the sample rate
    #	EMD_CONFIG controls how the EMD algorithm operations
    #         EMD_CONFIG$MAX_SIFT how many times the IMFs can be sifted
    #         EMD_CONFIG$MAX_IMF maximum number of IMFs that can be generated
    #         EMD_CONFIG$TOL Sifting stop criterion.
    #         EMD_CONFIG$STOP_RULE Make sure to read section on stop rules and make sure you understand what they imply!
    #         EMD_CONFIG$BOUNDARY How the start and stop of the time series are handled duing the spline process.
    #         EMD_CONFIG$SM Spline smoothing
    #         EMD_CONFIG$SPAR Smoothing parameter (only needed if sm is not none)
    #         EMD_CONFIG$WEIGHT Weight if "sm" is "spline"
    #         VERBOSE lets you know how many IMFs have been extracted.

    #OUTPUT is a list containing the original signal, IMFs, and information on EMD parameters.
    #Danny Bowman
    #UNC Chapel Hill
    
    #NOTES ON STOP RULES
    #These are rules to determine when sifting should stop and the result be saved as an IMF
    #The rules are implemented in the dcb_extractimf, which is taken almost verbatim from the extractimf function in the EMD package.
    #Here is how they work
    #TYPE1
    #type 1 finds the mean of the envelopes as defined by the spline function.  See Fig. 3b in Huang et al 1998:
    #The envelope is the dotted lines, the mean is the solid line.  In an IMF, the authors of the EMD function
    #think that this mean should be close to 0 (hence, below the defined tolerance).  This will force the IMF to be
    #symmetric, but I am concerned it may err towards frequency modulation at the expense of amplitude.
    #This criterion is not mentioned in literature I have read and they only cite Huang et al 1998 where
    #it is certainly not mentioned at all.  Instead they have derived this criterion from the definition of an IMF:
    #that at every point the mean defined by the local maximum and minimum envelope be zero.
    #TYPE2
    #Type 2 is the method used by Huang et al 1998.  It finds the standard deviation between two consecutive siftings;
    #when the SD falls below a certain tolerance value the sifting stops.  Huang et al recommend a tolerance value of 0.2 or 0.3
    #TYPE3
    #Type 3 is my addition to the EMD package.  It implements the S stoppage criterion per Huang et al 2008 (further described in
    #Huang et al 2003).  The S criterion means the sifting process stops when there has been no change in the number of extrema
    #and zero crossings for S iterations.  S is arbitrary but Huang's tests indicate around 3-8 works well but this should be tested
    #per Huang et al 2003 for each type of data otherwise it is somewhat arbitrary.

    tt=seq_len(length(sig))*dt
    tt=tt[which(!is.na(sig))]
    sig=sig[which(!is.na(sig))]
    emd_result=emd(sig, tt, max.sift=emd_config$max_sift, stoprule=emd_config$stop_rule, tol=emd_config$tol, 
        boundary=emd_config$boundary,sm=emd_config$sm,spar=emd_config$spar,weight=emd_config$weight, 
        check=FALSE, plot.imf=FALSE,max.imf=emd_config$max_imf)
    emd_result$original_signal=sig
    emd_result$dt=dt
    for(pname in names(emd_config))
    {
        emd_result[pname]=emd_config[pname]
    }
    invisible(emd_result)
}
