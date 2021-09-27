/**
 * File              : run.C
 * Author            : Anton Riedel <anton.riedel@tum.de>
 * Date              : 01.06.2021
 * Last Modified Date: 27.09.2021
 * Last Modified By  : Anton Riedel <anton.riedel@tum.de>
 */

Int_t run() {

  TFile *file = new TFile("ControlHistograms.root", "RECREATE");

  // pdf for multiplicity flucutations
  Int_t NumberOfParticles0 = 800;
  Int_t NumberOfParticles1 = 1200;
  TF1 *PDFmultiplicity =
      new TF1("Multiplicity", "1", NumberOfParticles0, NumberOfParticles1);
  PDFmultiplicity->Write();

  // pdf for phi distribution
  // set flow harmonics
  const Int_t N = 6;
  std::vector<Double_t> harmonics = {};
  for (int i = 0; i < N; ++i) {
    harmonics.push_back(0.05 + i * 0.01);
  }
  // write the formula
  // f(phi) = 1+2*Sum_i^N v_i cos(i*(phi-psi_i))
  TString formula = "1+2*(";

  for (int i = 0; i < N; i++) {
    formula += Form("[%d]*TMath::Cos(%d*(x-[%d]))", 2 * i, i + 1, 2 * i + 1);
    if (i != N - 1) {
      formula += TString("+");
    }
  }

  formula += ")";

  std::cout << formula << std::endl;

  Double_t phi0 = 0;
  Double_t phi1 = TMath::TwoPi();
  TF1 *PDFphi = new TF1("phi", formula, phi0, phi1);

  // set flow harmonics
  // symmetry planes will be set run by run
  for (std::size_t i = 0; i < harmonics.size(); i++) {
    PDFphi->SetParameter(2 * i, harmonics.at(i));
  }
  PDFphi->Write();

  // pdf for pt distribution
  Double_t pt0 = 0;
  Double_t pt1 = 5;
  TF1 *PDFpt =
      new TF1("pt", "x*TMath::Exp(-TMath::Sqrt(x^2+[0]^2)/[1])", pt0, pt1);
  Double_t m = 0.134;
  Double_t T = 0.44;
  PDFpt->SetParameter(0, m);
  PDFpt->SetParameter(1, T);
  PDFpt->Write();

  // pdf for eta distribution
  Double_t eta0 = -1.;
  Double_t eta1 = 1.;
  TF1 *PDFeta = new TF1("eta", "1", eta0, eta1);
  PDFeta->Write();

  // acceptance/weight histograms for phi
  Int_t binsPhi = 720;
  Double_t phi2 = TMath::TwoPi() / 6;
  Double_t phi3 = TMath::TwoPi() / 3;
  Double_t phiAccetance = 0.4;
  TH1D *AcceptancePhi =
      new TH1D("acceptance_phi", "acceptance_phi", binsPhi, phi0, phi1);
  TH1D *WeightPhi = dynamic_cast<TH1D *>(AcceptancePhi->Clone("weight_phi"));
  for (int i = 0; i < AcceptancePhi->GetNbinsX(); i++) {
    if (phi2 < AcceptancePhi->GetBinCenter(i + 1) &&
        phi3 > AcceptancePhi->GetBinCenter(i + 1)) {
      AcceptancePhi->SetBinContent(i + 1, phiAccetance);
      WeightPhi->SetBinContent(i + 1, 1. / phiAccetance);
    } else {
      AcceptancePhi->SetBinContent(i + 1, 1.);
      WeightPhi->SetBinContent(i + 1, 1.);
    }
  }
  AcceptancePhi->Write();
  WeightPhi->Write();

  // acceptance/weight histograms for pt
  Int_t binsPt = 1000;
  Double_t pt2 = 0.4;
  Double_t pt3 = 1.2;
  Double_t ptAccetance = 0.5;
  TH1D *AcceptancePt =
      new TH1D("acceptance_pt", "acceptance_pt", binsPt, pt0, pt1);
  TH1D *WeightPt = dynamic_cast<TH1D *>(AcceptancePt->Clone("weight_pt"));
  for (int i = 0; i < AcceptancePt->GetNbinsX(); i++) {
    if (pt2 < AcceptancePt->GetBinCenter(i + 1) &&
        pt3 > AcceptancePt->GetBinCenter(i + 1)) {
      AcceptancePt->SetBinContent(i + 1, ptAccetance);
      WeightPt->SetBinContent(i + 1, 1. / ptAccetance);
    } else {
      AcceptancePt->SetBinContent(i + 1, 1.);
      WeightPt->SetBinContent(i + 1, 1.);
    }
  }
  AcceptancePt->Write();
  WeightPt->Write();

  // acceptance/weight histograms for pt
  Int_t binsEta = 1000;
  Double_t eta2 = -0.1;
  Double_t eta3 = 0.3;
  Double_t etaAcceptance = 0.6;
  TH1D *AcceptanceEta =
      new TH1D("acceptance_eta", "acceptaance_eta", binsEta, eta0, eta1);
  TH1D *WeightEta = dynamic_cast<TH1D *>(AcceptanceEta->Clone("weight_eta"));
  for (int i = 0; i < AcceptanceEta->GetNbinsX(); i++) {
    if (eta2 < AcceptanceEta->GetBinCenter(i + 1) &&
        eta3 > AcceptanceEta->GetBinCenter(i + 1)) {
      AcceptanceEta->SetBinContent(i + 1, etaAcceptance);
      WeightEta->SetBinContent(i + 1, 1. / etaAcceptance);
    } else {
      AcceptanceEta->SetBinContent(i + 1, 1.);
      WeightEta->SetBinContent(i + 1, 1.);
    }
  }
  AcceptanceEta->Write();
  WeightEta->Write();

  // correlators we want to compute
  std::vector<std::vector<Int_t>> correlators = {{-2, 2},
                                                 {-5, -1, 6},
                                                 {-3, -2, 2, 3},
                                                 {-5, -4, 3, 3, 3},
                                                 {-2, -2, -1, -1, 3, 3},
                                                 {-6, -5, -1, 1, 2, 3, 6},
                                                 {-6, -6, -5, 2, 3, 3, 4, 5}};

  // compute theory values
  std::vector<Double_t> theoryValue;
  Double_t value;
  for (std::size_t i = 0; i < correlators.size(); i++) {
    value = 1.;
    for (std::size_t j = 0; j < correlators.at(i).size(); j++) {
      value *= harmonics.at(std::abs(correlators.at(i).at(j)) - 1);
    }
    theoryValue.push_back(value);
  }

  // custom seed
  UInt_t Seed = 2020201;

  // configure task for uniform acceptance
  AliAnalysisTaskAR *taskUA =
      new AliAnalysisTaskAR("UniformAcceptance", kFALSE);
  taskUA->SetMCOnTheFly(kTRUE);
  taskUA->SetUseCustomSeed(Seed);
  taskUA->SetCorrelators(correlators);
  taskUA->SetTrackControlHistogramBinning(kPT, binsPt, pt0, pt1);
  taskUA->SetTrackControlHistogramBinning(kPHI, binsPhi, phi0, phi1);
  taskUA->SetTrackControlHistogramBinning(kETA, binsEta, eta0, eta1);
  taskUA->SetEventControlHistogramBinning(
      kMUL, (NumberOfParticles1 - NumberOfParticles0), NumberOfParticles0,
      NumberOfParticles1);
  taskUA->SetEventControlHistogramBinning(kMULQ, NumberOfParticles1, 0,
                                          NumberOfParticles1);
  taskUA->SetEventControlHistogramBinning(kMULW, NumberOfParticles1, 0,
                                          NumberOfParticles1);

  // set pdfs
  taskUA->SetMCMultiplicityPdf(PDFmultiplicity);
  taskUA->SetMCKinematicPdf(kPT, PDFpt);
  taskUA->SetMCKinematicPdf(kETA, PDFeta);

  // configure task with non uniform acceptance and without weights
  AliAnalysisTaskAR *taskNUA =
      dynamic_cast<AliAnalysisTaskAR *>(taskUA->Clone("NonUniformAcceptance"));
  taskNUA->SetAcceptanceHistogram(kPHI, AcceptancePhi);
  taskNUA->SetAcceptanceHistogram(kPT, AcceptancePt);
  taskNUA->SetAcceptanceHistogram(kETA, AcceptanceEta);

  AliAnalysisTaskAR *taskNUAWW = dynamic_cast<AliAnalysisTaskAR *>(
      taskNUA->Clone("NonUniformAcceptanceWithWeights"));
  taskNUAWW->SetWeightHistogram(kPHI, WeightPhi);
  taskNUAWW->SetWeightHistogram(kPT, WeightPt);
  taskNUAWW->SetWeightHistogram(kETA, WeightEta);

  // run analysis
  taskUA->UserCreateOutputObjects();
  taskNUA->UserCreateOutputObjects();
  taskNUAWW->UserCreateOutputObjects();

  // loop over all events
  Int_t NumberOfEvents = 1000;
  Double_t Psi = 0.;
  for (int i = 0; i < NumberOfEvents; ++i) {

    std::cout << "Processing Event " << i << std::endl;

    // randomize psi
    Psi = gRandom->Uniform(TMath::TwoPi());
    for (std::size_t i = 0; i < harmonics.size(); i++) {
      PDFphi->SetParameter(2 * i + 1, Psi);
    }
    taskUA->SetMCKinematicPdf(kPHI, PDFphi);
    taskNUA->SetMCKinematicPdf(kPHI, PDFphi);
    taskNUAWW->SetMCKinematicPdf(kPHI, PDFphi);

    taskUA->UserExec(nullptr);
    taskNUA->UserExec(nullptr);
    taskNUAWW->UserExec(nullptr);
  }
  taskUA->Terminate(nullptr);
  taskNUA->Terminate(nullptr);
  taskNUAWW->Terminate(nullptr);

  taskNUAWW->GetEventControlHistogram(kSIM, kMUL, kAFTER)->Write();
  taskNUAWW->GetEventControlHistogram(kSIM, kMULQ, kAFTER)->Write();
  taskNUAWW->GetEventControlHistogram(kSIM, kMULW, kAFTER)->Write();
  taskNUAWW->GetTrackControlHistogram(kSIM, kPT, kBEFORE)->Write();
  taskNUAWW->GetTrackControlHistogram(kSIM, kPHI, kBEFORE)->Write();
  taskNUAWW->GetTrackControlHistogram(kSIM, kETA, kBEFORE)->Write();
  taskNUAWW->GetTrackControlHistogram(kSIM, kPT, kAFTER)->Write();
  taskNUAWW->GetTrackControlHistogram(kSIM, kPHI, kAFTER)->Write();
  taskNUAWW->GetTrackControlHistogram(kSIM, kETA, kAFTER)->Write();

  TList *outputList1 = taskUA->GetFinalResultProfilesList();
  TList *outputList2 = taskNUA->GetFinalResultProfilesList();
  TList *outputList3 = taskNUAWW->GetFinalResultProfilesList();
  // outputList3->Write();

  file->Close();

  std::vector<Double_t> corValueUA;
  std::vector<Double_t> corValueNUA;
  std::vector<Double_t> corValueNUAWW;

  TList *list;
  for (auto List : *outputList1) {
    list = dynamic_cast<TList *>(List);
    corValueUA.push_back(
        dynamic_cast<TProfile *>(list->At(0))->GetBinContent(1));
  }
  for (auto List : *outputList2) {
    list = dynamic_cast<TList *>(List);
    corValueNUA.push_back(
        dynamic_cast<TProfile *>(list->At(0))->GetBinContent(1));
  }
  for (auto List : *outputList3) {
    list = dynamic_cast<TList *>(List);
    corValueNUAWW.push_back(
        dynamic_cast<TProfile *>(list->At(0))->GetBinContent(1));
  }

  for (std::size_t i = 0; i < corValueUA.size(); i++) {
    std::cout << "####################" << std::endl;
    std::cout << "Theory:" << theoryValue.at(i) << std::endl;
    std::cout << "Uniform:" << corValueUA.at(i) << std::endl;
    std::cout << "Non-Uniform:" << corValueNUA.at(i) << std::endl;
    std::cout << "Non-Uniform WW:" << corValueNUAWW.at(i) << std::endl;
    std::cout << "####################" << std::endl;
  }

  // std::vector<Double_t> X(corValueUA.size());
  // std::iota(std::begin(X), std::end(X), 1);

  // // TGraph *theory = new TGraph(N, X, TheoryValues);
  // TGraph *MCUniform = new TGraph(Int_t(X.size()), X.data(),
  // corValueUA.data()); TGraph *MCNonUniform =
  //     new TGraph(Int_t(X.size()), X.data(), corValueNUA.data());

  // TH1F *hr = c->DrawFrame(0, -0.8, 7, 1);
  // hr->GetXaxis()->Set(N, 0, N);
  // for (int i = 0; i < N; ++i) {
  //   hr->GetXaxis()->SetBinLabel(i + 1, Form("#LT#LT%d#GT#GT", i + 2));
  // }
  // hr->GetYaxis()->SetTitle("#LT#LTk#GT#GT/10^{-k}");

  // MCUniform->Draw();
  // MCNonUniform->Draw("SAME");
  // c->SaveAs("result.pdf");

  return 0;
}
