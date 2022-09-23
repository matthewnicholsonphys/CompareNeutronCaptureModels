#include <string>
#include <vector>
#include <iostream>

#include "TTree.h"
#include "TVector3.h"
#include "TFile.h"
#include "THStack.h"
#include "TH1D.h"

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

  virtual std::vector<int> GetMultiplicity() const {
    std::vector<int> result;
    for (const auto& event : vec_of_events){result.push_back(event.vec_of_capt_cand.size());}
    return result;
  }

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
    
    const auto numb_of_events = sk2p2_tree->GetEntries();
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

    const auto numb_of_events = candidate_tree->GetEntries();
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
};

class ModelComparer {
private:
  std::vector<ModelRun*> models;
  std::string output_fname;
public:
  ModelComparer(const std::vector<ModelRun*>& m) : models{m} {}
  void Compare() const {

    THStack reco_mult_stack = CompareBetweenModels([](const ModelRun* m){return m->GetMultiplicity();},
						  "reco_mult_stack",
						  "Comparison Of Neutron Tagging Models - Reconstructed Multiplicity;Multiplicity;",
						  {11,0,10});
    reco_mult_stack.SaveAs("test_lambda.root");
						  

  }  

  template <typename ModelProperty>
  THStack CompareBetweenModels(const ModelProperty& thing_to_compare_of, const std::string&& stack_name, const std::string&& stack_title, const std::vector<double>&& bins) const {
    int colour_numb{0};
    THStack the_stack(stack_name.c_str(), stack_title.c_str());
    
    for (const auto& model : models){
      TH1D* h = new TH1D(model->GetModelName().c_str(), model->GetModelName().c_str(), bins.at(0), bins.at(1), bins.at(2));
      // root memory management sucks - THStack handles the deletion of this.
      h->SetFillColor(colour_numb+=2);
      h->SetLineWidth(5);
      for (const auto& c : thing_to_compare_of(model)){
	h->Fill(c);
      }
      the_stack.Add(h);
    }
    return the_stack;
  } 
};

int main(){

  const std::string ntag_file{"/raid6/marcus/ntag_cnn_skg4_ibd.root"};
  const std::string bdt_file{"/raid6/marcus/skg4ibd_bdtOut.root"};

  std::vector<ModelRun*> vec_of_models;

  NTag_ModelRun ntag_modelrun(ntag_file);
  vec_of_models.push_back(&ntag_modelrun);

  BDT_ModelRun bdt_modelrun(bdt_file);  
  vec_of_models.push_back(&bdt_modelrun);
 
  ModelComparer comparer(vec_of_models);
  comparer.Compare();
 
  return 0;
}
