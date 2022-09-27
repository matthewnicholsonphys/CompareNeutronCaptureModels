#include <string>
#include <vector>
#include <iostream>

#include "TTree.h"
#include "TVector3.h"
#include "TFile.h"
#include "THStack.h"
#include "TH1D.h"

#include "EventTrueCaptures.h"

struct CaptureCandidate {
  double likelihood_metric{0};
  int matched{0};  // 0=mistag, 1=H, 2=Gd, 3=decay e                            
  double reco_t{0};
  double true_t{0};
  double t_err{0};
  TVector3 reco_vtx;
  TVector3 true_vtx;
  double vtx_err{0};
  double dwall{0};
  double true_gamma_energy{0};
  int eventnum{0};
};
  
struct Event {
  int event_number;
  std::vector<CaptureCandidate> vec_of_capt_cand;
  
};

class ModelRun {
private:
protected:
  std::string model_name;
  std::string data_fname;
  std::vector<Event> vec_of_events;
public:
  ModelRun(const std::string& f) : data_fname{f} {
    std::cout << "ModelRun constructor starts\n";

    /*

    TFile* data_file = TFile::Open(data_fname.c_str());
    TTree* mc_tree = (TTree*) data_file->Get("mc");
    
    auto numb_of_entries = mc_tree->GetEntries();
    for (auto entry_idx = 0; entry_idx < numb_of_entries; ++entry_idx){
      mc_tree->GetEntry(entry_idx);
    }
    data_file->Close();

    */

    std::cout << "ModelRun constructor ends\n";
  }

  virtual std::vector<int> GetNumbOfCandidates() const {
    std::vector<int> result;
    for (const auto& event : vec_of_events){result.push_back(event.vec_of_capt_cand.size());}
    return result;
  }

  virtual std::vector<double> GetLikelihoodMetrics() const {
    std::vector<double> result;
    for (auto&& event : vec_of_events){
      for (auto&& cap_cand : event.vec_of_capt_cand){
	result.push_back(cap_cand.likelihood_metric);
      }
    }
    return result;
  }

  virtual std::vector<int> GetReconstructedMultiplicity() const = 0;

  const std::string& GetModelName() const {return model_name;}
};

class BDT_ModelRun : public ModelRun {
private:

public:
  BDT_ModelRun(const std::string& f) : ModelRun{f} {
    std::cout << "BDT_ModelRun constructor starts\n";
    model_name = "BDT";
    TFile* data_file = TFile::Open(data_fname.c_str());
    TTree* sk2p2_tree = (TTree*) data_file->Get("sk2p2");
    int numb_of_candidates;
    float candidate_metrics[500];
    float candidate_xs[500];
    float candidate_ys[500];
    float candidate_zs[500];
    float candidate_ts[500];
    sk2p2_tree->SetBranchAddress("np",&numb_of_candidates);
    sk2p2_tree->SetBranchAddress("neutron5",&candidate_metrics);  // 1 = neutron, 0 = not neutron (wiki)                                
    // check all the following, not sure of variable names                                                                           
    sk2p2_tree->SetBranchAddress("nvx",&candidate_xs); // total guesses that these are the branches we need.                            
    sk2p2_tree->SetBranchAddress("nvy",&candidate_ys);
    sk2p2_tree->SetBranchAddress("nvz",&candidate_zs);
    sk2p2_tree->SetBranchAddress("dt",&candidate_ts);  // dt positron->neutron                                                          
    // these are offset of 20E3 and need scaling by /1E3? positron time should be 0 for neutron gun...
    
    const auto numb_of_events = 48667; //sk2p2_tree->GetEntries();
    for (auto event_idx = 0; event_idx < numb_of_events; ++event_idx){
      sk2p2_tree->GetEntry(event_idx);
      Event event;
      event.event_number = event_idx;
      for (auto cand_idx = 0; cand_idx < numb_of_candidates; ++cand_idx){
	CaptureCandidate cap_cand;
	cap_cand.likelihood_metric = candidate_metrics[cand_idx];
	cap_cand.reco_t = candidate_ts[cand_idx];
	cap_cand.reco_vtx = TVector3(candidate_xs[cand_idx],
				     candidate_ys[cand_idx],
				     candidate_zs[cand_idx]);
	cap_cand.matched = 0;

	event.vec_of_capt_cand.push_back(cap_cand);
      }
      vec_of_events.push_back(event);
    }
    data_file->Close();
    std::cout << "BDT_ModelRun constructor ends\n";
  }

  std::vector<int> GetReconstructedMultiplicity() const override {
    std::vector<int> result;
    for (auto&& event : vec_of_events){
      int numb_of_neutrons{0};
      for(auto&& cap_cand : event.vec_of_capt_cand){
	if (cap_cand.likelihood_metric == 1){++numb_of_neutrons;}
      }
      result.push_back(numb_of_neutrons);
    }
    return result;
  } 

};

class NTag_ModelRun : public ModelRun {
private:

public:
  NTag_ModelRun(const std::string& f) : ModelRun{f} {
    std::cout << "NTag_ModelRun constructor starts\n";
    model_name = "NTAG";
    TFile* data_file = TFile::Open(data_fname.c_str());
    TTree* candidate_tree = (TTree*) data_file->Get("candidates");
    std::vector<float> candidate_times;
    std::vector<float> candidate_metrics;
    std::vector<float> truths;
    std::vector<float> candidate_xs;
    std::vector<float> candidate_ys;
    std::vector<float> candidate_zs;

    std::vector<float>* candidate_times_ptr = &candidate_times;
    std::vector<float>* candidate_metrics_ptr = &candidate_metrics;
    std::vector<float>* truths_ptr = &truths;
    std::vector<float>* candidate_xs_ptr = &candidate_xs;
    std::vector<float>* candidate_ys_ptr = &candidate_ys;
    std::vector<float>* candidate_zs_ptr = &candidate_zs;

    candidate_tree->SetBranchAddress("ReconCT",&candidate_times_ptr);
    // ReconCT=mean, TRMS=RMS of candidate hit Ts. Note both need to be scaled by 1000                                               
    // to compare to a TrueCapture time...??                                                                                         
    candidate_tree->SetBranchAddress("TMVAOutput",&candidate_metrics_ptr);
    // 'captureType' is capture nucleus: 0=Mistag, 1=H, 2=Gd, 3=decay e-                                                             
    candidate_tree->SetBranchAddress("CaptureType",&truths_ptr);
    //not usually saved to file:
    candidate_tree->SetBranchAddress("TrmsFitVertex_X",&candidate_xs_ptr);
    candidate_tree->SetBranchAddress("TrmsFitVertex_Y",&candidate_ys_ptr);
    candidate_tree->SetBranchAddress("TrmsFitVertex_Z",&candidate_zs_ptr);

    const auto numb_of_events = 48667; //candidate_tree->GetEntries();
    for (auto event_idx = 0; event_idx < numb_of_events; ++event_idx){
      candidate_tree->GetEntry(event_idx);
      Event event;
      event.event_number = event_idx;
      const auto numb_of_candidates_in_event = candidate_times.size();
      for (auto cand_idx = 0; cand_idx < numb_of_candidates_in_event; ++cand_idx){
	CaptureCandidate cap_cand;

	cap_cand.likelihood_metric = candidate_metrics.at(cand_idx);
	cap_cand.reco_t = candidate_times.at(cand_idx);
	cap_cand.reco_vtx = TVector3(candidate_xs.at(cand_idx), 
				     candidate_ys.at(cand_idx),
				     candidate_zs.at(cand_idx));
	cap_cand.matched = truths.at(cand_idx);

	event.vec_of_capt_cand.push_back(cap_cand);
      }
      vec_of_events.push_back(event);
    }
    data_file->Close();
    std::cout << "NTag_ModelRun constructor ends\n";
  }

  std::vector<int> GetReconstructedMultiplicity() const override {
    std::vector<int> result;
    for (auto&& event : vec_of_events){
      int numb_of_neutrons{0};
      for(auto&& cap_cand : event.vec_of_capt_cand){
	if (cap_cand.likelihood_metric != 0){++numb_of_neutrons;}
      }
      result.push_back(numb_of_neutrons);
    }
    return result;
  } 
  
};

class ModelComparer {
private:
  std::vector<ModelRun*> models;

  template <typename ModelPropLmbda1, typename ModelPropLmbda2>
  THStack CompareProperties(const ModelRun* model_run,
			    const ModelPropLmbda1 prop1_of, std::string&& prop1_name,
			    const ModelPropLmbda2 prop2_of, std::string&& prop2_name,
			    std::string&& stack_name,
			    std::string&& stack_title,
			    std::array<int, 3>&& bins){
    int colour{0};

    THStack the_stack(stack_name.c_str(), stack_title.c_str());
    
    TH1D* hist1 = new TH1D((prop1_name + "_" + model_run->GetModelName()).c_str(),
			   prop1_name.c_str(),
			   bins.at(0), bins.at(1), bins.at(2));
    
    hist1->SetFillColor(colour+=2);
    hist1->SetLineWidth(4);
    for (auto&& elem : prop1_of(model_run)){
      hist1->Fill(elem);
    }
    the_stack.Add(hist1);

    TH1D* hist2 = new TH1D((prop2_name + "_" + model_run->GetModelName()).c_str(),
			   prop2_name.c_str(),
			   bins.at(0), bins.at(1), bins.at(2));
    
    hist2->SetFillColor(colour+=2);
    hist2->SetLineWidth(4);
    for (auto&& elem : prop2_of(model_run)){
      hist2->Fill(elem);
    }
    the_stack.Add(hist2);

    return the_stack;
  }

  template <typename ModelPropertyLmbda>
  THStack CompareBetweenModels(const ModelPropertyLmbda& return_property_of,
			       std::string&& stack_name,
			       std::string&& stack_title,
			       std::array<int, 3>&& bins){
  int colour{0};

  THStack the_stack(stack_name.c_str(), stack_title.c_str());
  for (const auto& model : models){
    //root memory managment sucks - THStack handles deletion.
    TH1D* hist = new TH1D((stack_name + "_" + model->GetModelName()).c_str(), 
			  model->GetModelName().c_str(), 
			  bins.at(0), bins.at(1), bins.at(2));
    hist->SetFillColor(colour+=2);
    hist->SetLineWidth(5);
    for (auto&& elem : return_property_of(model)){
      hist->Fill(elem);
    }
    the_stack.Add(hist);
  }
  return the_stack;
}

public:
  ModelComparer(std::vector<ModelRun*> m) : models{m} {}

  void CompareAndSaveToFile(std::string&& out_fname){
    TFile* fout = new TFile(out_fname.c_str(), "recreate");

    THStack numb_of_cand_stack = CompareBetweenModels([](const ModelRun* m){return m->GetNumbOfCandidates();},
						       "numb_of_cand_stack",
						       "Comparison Of Neutron Models - Number Of Candidates; numb of candidates",
						       {11, 0, 10});
    numb_of_cand_stack.Write();

    THStack likelihood_metrics_stack = CompareBetweenModels([](const ModelRun* m){return m->GetLikelihoodMetrics();},
						      "likelihood_metrics_stack",
						      "Comparison Of Neutron Models - Likelihood Metrics; metric",
						      {11, 0, 10});
    likelihood_metrics_stack.Write();


    THStack reco_mult_stack = CompareBetweenModels([](const ModelRun* m){return m->GetReconstructedMultiplicity();},
						   "reco_mult_stack",
						   "Comparison Of Neutron Models - Reconstructed Multiplicity; multiplicity",
						   {11, 0, 10});
    reco_mult_stack.Write();
						   



    fout->Close();
    delete fout;

  }
};

int main(){

  TClass::GetClass("TVector3")->IgnoreTObjectStreamer();

  const std::string ntag_file{"/raid6/marcus/ntag_cnn_skg4_ibd.root"};
  const std::string bdt_file{"/raid6/marcus/skg4ibd_bdtOut.root"};

  std::vector<ModelRun*> vec_of_models;

  NTag_ModelRun ntag_modelrun(ntag_file);
  vec_of_models.push_back(&ntag_modelrun);

  BDT_ModelRun bdt_modelrun(bdt_file);  
  vec_of_models.push_back(&bdt_modelrun);
 
  ModelComparer comparer(vec_of_models);
  comparer.CompareAndSaveToFile("./neutron_tagging_comparison_output.root");
 
  return 0;
}
