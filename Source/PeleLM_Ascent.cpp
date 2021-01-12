//
// Demonstrates how to generate of a Conduit Mesh Blueprint 
// description of an AMReX Single Level dataset and render this
// data in situ using with ALPINE Ascent.
// 

#include <AMReX_Utility.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Conduit_Blueprint.H>

#include <PeleLM.H>

#include <conduit/conduit.hpp>
#include <conduit/conduit_blueprint.hpp>
#include <conduit/conduit_relay.hpp>

#include <ascent.hpp>

#include "AMReX_buildInfo.H"

using std::string;
using namespace amrex;

using namespace conduit;
using namespace ascent;

void 
PeleLM::goAscent (int nstep)
{
  amrex::Print() << " goAscent nstep " << nstep << std::endl;
  std::vector<std::pair<int,int> > plot_var_map;
  Vector<std::string> var_names;
  for (int typ = 0; typ < desc_lst.size(); typ++)
  {
     for (int comp = 0; comp < desc_lst[typ].nComp();comp++)
     {
        if (parent->isStatePlotVar(desc_lst[typ].name(comp)) && 
            desc_lst[typ].getType() == IndexType::TheCellType())
        {
          plot_var_map.push_back(std::pair<int,int>(typ,comp));
          var_names.push_back(desc_lst[typ].name(comp));
          amrex::Print()<<"Var "<<desc_lst[typ].name(comp)<<"\n";
        }
     }
  }

  int n_data_items = plot_var_map.size();
  amrex::Print() << " n_data_items " << n_data_items << std::endl;

  Real time = state[State_Type].curTime();

  // We combine all of the multifabs -- state, derived, etc -- into one
  // multifab -- plotMF.
  // NOTE: we are assuming that each state variable has one component,
  // but a derived variable is allowed to have multiple components.
        amrex::Print()<<"0000000000000000000000000000000000000000000000000000000000000000000000\n";
        //const std::list<std::string>& plot_vars = amrptr->statePlotVars();
        Vector<const MultiFab*> mfs;
        Vector<Geometry> geoms;
        const int nGrow = 1;
        MultiFab  plotMF(grids,dmap,1,nGrow,MFInfo(),Factory());
        Vector<int> level_steps;
        Vector<IntVect> ref_ratios;

        for(int lev = 0; lev <= parent->finestLevel(); ++lev)
        {
          MultiFab*  level_mf = new MultiFab(parent->getLevel(lev).boxArray(),
                                             parent->getLevel(lev).DistributionMap(),
                                             n_data_items,
                                             nGrow,
                                             MFInfo(),
                                             parent->getLevel(lev).Factory());
          level_mf->FillBoundary(); // Added by Matt for rendering
          mfs.push_back(level_mf);
          MultiFab* this_dat = 0;
          //
          // Cull data from state variables --
          //
          int cnt = 0;
          for (int i = 0; i < static_cast<int>(plot_var_map.size()); i++)
          {
            int type  = plot_var_map[i].first;
            int comp = plot_var_map[i].second;
            this_dat = &parent->getLevel(lev).get_new_data(type);
            MultiFab::Copy(*level_mf,*this_dat,comp,i,1,nGrow);
            //amrex::Print()<< "Copy cnt " << type << " " << comp << std::endl;
            cnt++;
          }

          geoms.push_back(parent->getLevel(lev).Geom());
          level_steps.push_back(parent->levelSteps(lev));
          ref_ratios.push_back(parent->getLevel(lev).fineRatio());

        }

  //MultiFab  plotMF(grids,dmap,1,nGrow,MFInfo(),Factory());
  //MultiFab  plotMF(grids,dmap,n_data_items,nGrow,MFInfo(),Factory());

  /////////////////////////////
  // Setup Ascent
  /////////////////////////////
  // Create an instance of Ascent
  Ascent ascent;
  Node open_opts;

  // for the MPI case, provide the mpi comm
#ifdef BL_USE_MPI
  open_opts["mpi_comm"] = MPI_Comm_c2f(ParallelDescriptor::Communicator());
#endif

  ascent.open(open_opts); 
  ///////////////////////////////////////////////////////////////////
  // Wrap our AMReX Mesh into a Conduit Mesh Blueprint Tree
  ///////////////////////////////////////////////////////////////////

  // Write a plotfile of the current data and write a conduit blueprint
  // file as well

  // in Base/AMReX_PlotFileUtil.cpp

  ///////////////////////////////////////////////////////////////////
  // Wrap our AMReX Mesh into a Conduit Mesh Blueprint Tree
  ///////////////////////////////////////////////////////////////////

        conduit::Node bp_mesh;
        MultiLevelToBlueprint( parent->finestLevel()+1,
                               mfs,
                               var_names,
                               geoms,
                               parent->cumTime(),
                               level_steps,
                               ref_ratios,
                               bp_mesh);

        conduit::Node verify_info;
        if(!conduit::blueprint::mesh::verify(bp_mesh,verify_info))
        {
          // verify failed, print error message
          ASCENT_INFO("Error: Mesh Blueprint Verify Failed!");
          // show details of what went awry
          verify_info.print();
        }
        else
        {
          amrex::Print()<< " everything A-ok" << std::endl;
        //      verify_info.print();
        }

  ///////////////////////////////////////////////////////////////////
  // Render with Ascent
  ///////////////////////////////////////////////////////////////////

  // add a scene with a pseudocolor plot
  Node scenes;
  scenes["s1/plots/p1/type"] = "pseudocolor";
  scenes["s1/plots/p1/field"] = "temperature";

  //Set the output file name (ascent will add ".png")
  const std::string& png_out = amrex::Concatenate("ascent_render_",n_data_items,5);
  scenes["s1/image_prefix"] = png_out;

  ///////////////////////////////////////////////////////////////////

  // setup actions
  Node actions;
  Node &add_act = actions.append();
  add_act["action"] = "add_scenes";
  add_act["scenes"] = scenes;

  // add_act["action"] = "add_extracts";
  // add_act["extracts"] = extracts;

  actions.append()["action"] = "execute";
  actions.append()["action"] = "reset";

#if 0
  conduit::Node add_pipelines = actions.append();
  add_pipelines["action"] = "add_pipelines";
  add_pipelines["pipelines"] = pipelines;
#endif

  ascent.publish(bp_mesh);

  ascent.execute(actions);

//  ascent.close();
}
