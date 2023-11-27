//Digitizer data from the LAB

struct slimport_data_t {
	ULong64_t	timetag; //time stamp
	UInt_t		baseline;
	UShort_t	qshort; //integration with shorter time
	UShort_t	qlong; //integration with longer time
	UShort_t	pur;
	UShort_t	samples[4096];
};

TH1D* getHistoFromTree(const char *name_file, int numBins, double minX, double maxX) {
	// variables
	slimport_data_t indata;
	TFile *infile = new TFile(name_file);
	TTree *intree = (TTree*)infile->Get("acq_tree_0");
	TBranch *inbranch = intree->GetBranch("acq_ch0");
	inbranch->SetAddress(&indata.timetag);
	TH1D *h_spectrum = new TH1D("h_spectrum","Total spectrum",numBins,minX,maxX);
	// histogram filling
	for (int i=0; i<inbranch->GetEntries(); i++) {
		inbranch->GetEntry(i);
		h_spectrum->Fill(indata.qlong);
	}
	// return
	return h_spectrum;
}

TH1D* getHistoForChannelFromTree(const char *name_file, short chan, int numBins, double minX, double maxX) {
	// variables
	slimport_data_t indata;
	TFile *infile = new TFile(name_file);
	TTree *intree = (TTree*)infile->Get("acq_tree_0");
	TBranch *inbranch = intree->GetBranch(Form("acq_ch%d",chan));
	inbranch->SetAddress(&indata.timetag);
	TH1D *h_spectrum = new TH1D("h_spectrum","Total spectrum",numBins,minX,maxX);
	// histogram filling
	for (int i=0; i<inbranch->GetEntries(); i++) {
		inbranch->GetEntry(i);
		h_spectrum->Fill(indata.qlong);
	}
	// return
	return h_spectrum;
}


TH1D* getHistoWithFilter(const char *name_file, int numBins, double minX, double maxX, double lowThr = 0, double highThr = 999999) {
	// variables
	slimport_data_t indata;
	TFile *infile = new TFile(name_file);
	TTree *intree = (TTree*)infile->Get("acq_tree_0");
	TBranch *inbranch = intree->GetBranch("acq_ch0");
	inbranch->SetAddress(&indata.timetag);
	TH1D *h_spectrum = new TH1D("h_spectrum","Total spectrum",numBins,minX,maxX);
	// histogram filling
	for (int i=0; i<inbranch->GetEntries(); i++) {
		inbranch->GetEntry(i);
		if (indata.qlong>lowThr && indata.qlong<highThr) {
			h_spectrum->Fill(indata.qlong);
		}
	}
	// return
	return h_spectrum;
}


TGraph *getSignal(const char *name_file, int nSamples=250, short chan=0, int nrEv=1){
	// variables
	slimport_data_t indata;
	TFile *infile = new TFile(name_file);
	TTree *intree = (TTree*)infile->Get("acq_tree_0");
	TBranch *inbranch = intree->GetBranch(Form("acq_ch%d",chan));
	inbranch->SetAddress(&indata.timetag);
	TGraph *graph = new TGraph();
	
	//Setting the desired event
	inbranch->GetEntry(nrEv);

	//Looping over the samples
	for (int i=0; i<nSamples; ++i){
		graph->SetPoint(i, i, indata.samples[i]);
	}
	return graph;
}

/////////////////////////////////////////////////////////////////////////////////////
//			ADDED PARTS AS OF 29/10,	~GIOELE


TBranch * getBranchFromTree(const char *name_file, const char *name_tree, const char *name_branch) {
	// variables
	slimport_data_t indata;
	TFile *infile = new TFile(name_file);
	TTree *intree = (TTree*)infile->Get(name_tree);
	TBranch *inbranch = intree->GetBranch(name_branch);

	return inbranch;
}

//if in doubt on what params use, leave default
TH1D* createSpectrum(TBranch* inbranch, bool ch1, int numBins = 1e3, double minX = 1, double maxX = -1, int nSamples = 250) {
	// variables
	slimport_data_t indata;
	//TFile *infile = new TFile(name_file);
	//TTree *intree = (TTree*)infile->Get("acq_tree_0");
	//TBranch *inbranch = intree->GetBranch("acq_ch0");
	inbranch->SetAddress(&indata.timetag);

    TH1D * spectrum = new TH1D("quantized_spectrum","Total spectrum",numBins,minX,maxX);
	//TH1D * spectrum = new TH1D();

    // cycle over all events, and then calculate integral
    for (int event = 0; event < inbranch -> GetEntries(); event++){
        inbranch -> GetEntry(event);

        //integration via a trapezoid, uniformely-discrete grid approximation
        double integral = 0;

		double baseline=0;
			for (int j = 20; j < 50 ; j++)
				baseline += indata.samples[j];
			baseline /= 30.0;

        // cycle over quantized samples of waveform and add area of infinitesimal trapezoid
        for(int i = 0; i < nSamples - 1; i++){
            //UShort_t f1,f2;
            //UShort_t dx = 1;    //grid spacing, CHECK IF RIGHT VALUE (1 channel)
			
			double f1,f2;
			double dx = 1;

            // access quantized values
            f1 = indata.samples[i];
            f2 = indata.samples[i+1];

            // calculate integral via trapezoid approx.
            integral += baseline - (f1 + f2) * dx / 2.0;
        }

		//convert from channels to energy
			// Double_t slope_D1 = 0.08148;
			// Double_t bias_D1 = -14.57;

			// Double_t slope_D2 = 0.08042;
			// Double_t bias_D2 = -28.52;

			
			Double_t slope_D1 = 0.08148;
			Double_t bias_D1 = -14.57;

			Double_t slope_D2 = 0.0798906;
			Double_t bias_D2 = -3.00144e1;
			
			if(ch1){
				integral = integral * slope_D1 - bias_D1;
			}else{
				integral = integral * slope_D2 - bias_D2;
			}
			


        // check thresholds and populate histogram
        //if (integral >= lowThr && integral <= highThr){
            spectrum -> Fill(integral);
            //if(integral > max_val) max_val = integral;
            //if(integral < min_val) min_val = integral;
        //}
    }

	// return
	return spectrum;
}

/////////////////////////////////////////////////////////////////////////////////////
//			ADDED PARTS AS OF 1/11,	~GIOELE
//			EDITED ON 4/11

TH1D* getTimeSpectrum(TBranch* ch0, TBranch* ch1, int D, float F, int nSamples=230, int initial_offset = 20, int numBins = 200, double l_thresh_energy=-1)	{
	TH1D * timeSpectrum = new TH1D("digitizer_t_spec", "Time spectrum from Digitizer", numBins, 1, 0);

	for ( int event = 0; event < ch0 -> GetEntries(); event++ ) {
		//		------------------ PART 0: LOAD DATA AND CFD ------------------
		slimport_data_t indata_ch0;
		slimport_data_t indata_ch1;

		//loading channels
		ch0 -> SetAddress(&indata_ch0.timetag);
		ch0 -> GetEvent(event);
		ch1 -> SetAddress(&indata_ch1.timetag);
		ch1 -> GetEvent(event);


		if (l_thresh_energy >= 0){
			//	------------------ PART 0.5: ENERGY THRESHOLD ------------------
			//convert from energy to channels			
			Double_t slope_D1 = 0.08148;
			Double_t bias_D1 = -14.57;

			Double_t slope_D2 = 0.08042;
			Double_t bias_D2 = -28.52;


			Double_t l_thresh_ch0 = (l_thresh_energy - bias_D1) / slope_D1;
			Double_t l_thresh_ch1 = (l_thresh_energy - bias_D2) / slope_D2;


			//integration via a trapezoid, uniformely-discrete grid approximation
			double integral_ch0 = 0;
			double integral_ch1 = 0;




			//		------------------ PART 1: ZERO CROSSING POINT ------------------
			//basic idea, first 50 bins of waveform are just noise around the baseline, take that and average out the noise
			double baseline_ch0 = 0;
			double baseline_ch1 = 0;
			for (int i = initial_offset; i < 50 ; i++) {
				baseline_ch0 += indata_ch0.samples[i];
				baseline_ch1 += indata_ch1.samples[i];
			}
			baseline_ch0 /= 50.0 - initial_offset;
			baseline_ch1 /= 50.0 - initial_offset;

			// for (int i = 0; i< nSamples; i++) {
			// 	indata_ch0.samples[i] -= baseline_ch0;
			// 	indata_ch1.samples[i] -= baseline_ch1;
			// }


			// cycle over quantized samples of waveform and add area of infinitesimal trapezoid
			for(int i = 0; i < nSamples - 1; i++){
				//UShort_t f1,f2;
				//UShort_t dx = 1;    //grid spacing, CHECK IF RIGHT VALUE (1 channel)
				
				double f1_ch0,f2_ch0;
				double f1_ch1,f2_ch1;
				double dx = 1.0;

				// access quantized values
				f1_ch0 = indata_ch0.samples[i];
				f2_ch0 = indata_ch0.samples[i+1];
				f1_ch1 = indata_ch1.samples[i];
				f2_ch1 = indata_ch1.samples[i+1];

				// calculate integral via trapezoid approx.
				integral_ch0 += baseline_ch0 * dx - (f1_ch0 + f2_ch0) * dx / 2.0;
				integral_ch1 += baseline_ch1 * dx - (f1_ch1 + f2_ch1) * dx / 2.0;
			}
			if( integral_ch0 < l_thresh_ch0 || integral_ch1 < l_thresh_ch1) continue;	//integral is negative due to electr. charge
		}



		//storage for the two channels' binomuials
		// TH1D * wave_ch0 = new TH1D("waveform_ch0", "waveform for ch0", nSamples, 1, 0); //putting xmin > xmax lets root choose the optimal range
		// TH1D * wave_ch1 = new TH1D("waveform_ch1", "waveform for ch1", nSamples, 1, 0);
		double binomial_ch0[nSamples - initial_offset];
		double binomial_ch1[nSamples - initial_offset];


		for ( int i = initial_offset - 1; i < nSamples; i++) {
			//CFD implementation
			binomial_ch0[i] = - indata_ch0.samples[i - D] + F * indata_ch0.samples[i];
			binomial_ch1[i] = - indata_ch1.samples[i - D] + F * indata_ch1.samples[i];
		}

		//		------------------ PART 1: ZERO CROSSING POINT ------------------
		//basic idea, first 50 bins of waveform are just noise around the baseline, take that and average out the noise
		double baseline_ch0 = 0;
		double baseline_ch1 = 0;
		for (int i = initial_offset; i < 50 ; i++) {
			baseline_ch0 += binomial_ch0[i];
			baseline_ch1 += binomial_ch1[i];
		}
		baseline_ch0 /= 50.0 - initial_offset;
		baseline_ch1 /= 50.0 - initial_offset;

		for (int i = 0; i< nSamples; i++) {
			binomial_ch0[i] -= baseline_ch0;
			binomial_ch1[i] -= baseline_ch1;
		}

		/////////////////////////////////////////////////////////////////
		// debugging: show plot of binomial signal
		// TGraph * draw_binomial_ch0 = new TGraph();
		// TGraph * draw_binomial_ch1 = new TGraph();

		// for ( int i = initial_offset; i < nSamples; ++i ) {
		// 	draw_binomial_ch0 -> SetPoint(i, i, binomial_ch0[i]);
		// 	draw_binomial_ch1 -> SetPoint(i, i, binomial_ch1[i]);
		// }

		// draw_binomial_ch0 -> Draw();
		// draw_binomial_ch1 -> Draw("same");
		/////////////////////////////////////////////////////////////////

		// finding the zero crossing point: scroll through the bipolar and find the point that goes from neg to pos
		double zcp_ch0 = -1;
		double zcp_ch1 = 1;
		for (int i = 0; i < 100; i++) { // limit the research to 100 since the waveform is in range (70, 150)

			// the zcp is such that the binomial is negative before it, positive after it and the rest of the signal is rapidly increasing
			if ( binomial_ch0[i] <= 0 && binomial_ch0[i+1] > 0 && binomial_ch0[i + 5] > binomial_ch0[i] + 5){	//the noise is in the range (-2,2)
				if ( binomial_ch0[i] == 0 ) zcp_ch0 = i;
				// not enough precision, switching to an interpolation
				//else zcp_ch0 = (double) i + 0.5;
				else zcp_ch0 = (double) i - binomial_ch0[i] / ( binomial_ch0[i + 1] - binomial_ch0[i] );
				break;	//stop searching
			}
		}

		for (int i = 0; i < 100; i++) { // limit the research to 100 since the waveform is in range (70, 150)
			//same thing for ch1
			if ( binomial_ch1[i] <= 0 && binomial_ch1[i+1] > 0 && binomial_ch1[i + 5] > binomial_ch1[i] + 5){
				if ( binomial_ch1[i] == 0 ) zcp_ch1 = i;
				//else zcp_ch1 = (double) i + 0.5;
				else zcp_ch1 = (double) i - binomial_ch1[i] / ( binomial_ch1[i + 1] - binomial_ch1[i] );
				break;	//stop searching
			}
		}
		double t_difference = zcp_ch0 - zcp_ch1; 
		timeSpectrum -> Fill(t_difference);
	}
	
	return timeSpectrum;
}


void setupGausnWLinBKG(){
	// fitting function
	auto gausWithLin = new TF1("gaus with linear bkg","gaus(0) + pol1(3)");
	gausWithLin->SetParLimits(0, 150, 30e3);
	gausWithLin->SetParLimits(1, 12, 14);
	gausWithLin->SetParLimits(2, 0.05, 2);
	gausWithLin->SetParLimits(3, 0, 1e3);
}

///////////////////////////////////////////////////////////////////////
//				ADDED PARTS AS OF 21/11		~ GIOELE
TH1D* GetCoincidence(const char *name_file, short chan, int numBins, double minX, double maxX) {
	// variables
	slimport_data_t indata;
	slimport_data_t trackers_sig;
	
	TFile *infile = new TFile(name_file);
	TTree *intree = (TTree*)infile->Get("acq_tree_0");
	
	TBranch *inbranch = intree->GetBranch(Form("acq_ch%d",chan));
	TBranch *trackers = intree->GetBranch("acq_ch3");

	inbranch->SetAddress(&indata);
	trackers->SetAddress(&trackers_sig);

	TH1D *h_spectrum = new TH1D("h_spectrum","Total spectrum",numBins,minX,maxX);
	// histogram filling

	Int_t trackers_entries = trackers->GetEntries();
	Int_t inbranch_entries = inbranch->GetEntries();
	
	for (int i = 0; i < trackers_entries; i++) {
		// inbranch->GetEntry(i);
		// trackers->GetEntry(i);
		
		// if ( trackers_sig.qlong > 1e3 && trackers_sig.qlong < 30e3 ){  
		// 	h_spectrum->Fill( indata.qlong );
		// 	//std::cout << "Trackers timetag: " << trackers_sig.timetag << "\n" << "Indata timetag: " << indata.timetag << endl;
		// 	//trackers->GetEntry(i+1);
		// 	//std::cout << "Successive Tracker qlong: " << trackers_sig.qlong << endl;
			
		// 	tdiff[counter] = trackers_sig.timetag - indata.timetag;
		// 	counter++;
		// }
		//	OLDER VERSION ^^^^^

		trackers->GetEntry(i);

		if ( trackers_sig.qlong > 1e3 && trackers_sig.qlong < 25e3 ){
			for ( int j=0; j < inbranch_entries; j++ ){
				// check timetags
				// if ( i < 1500 ) trackers->GetEntry(j);
				// else if ( i >= 1500 ) trackers->GetEntry(i-1500+j);
				inbranch->GetEntry(j);
				if ( trackers_sig.timetag == indata.timetag ) {
					h_spectrum->Fill( indata.qlong );
					break;
				}

				if ( j == inbranch_entries - 1 ) std::cout << "UNPAIRED EVENT!\n";
			}
		}
	}

	// return
	return h_spectrum;
}