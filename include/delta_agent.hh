//
// Interface of an simple agent.
//
// by Anton Crombach, A.B.M.Crombach at uu.nl
//

#ifndef _BLOWHOLE_DELTAAGENT_H_
#define _BLOWHOLE_DELTAAGENT_H_

#include "defs.hh"
#include "genreg_agent.hh"
#include "genome.hh"
#include "hop_graph.hh"

namespace blowhole {

    /// \class DeltaAgent
    /// \brief Agent with genome and network.
    ///
    /// The \c DeltaAgent has a fitness function that rewards being different
    /// from what you were in the previous environment.
    class DeltaAgent : public GenRegAgent {
        public:
        /// Default ctor.
        DeltaAgent();
        /// Ctor w/ args \c tag, \c genome, \c reference state
        DeltaAgent( uint, Genome *, const boost::dynamic_bitset<> & );
        /// Copy ctor.
        DeltaAgent( const DeltaAgent & );
        /// Destructor.
        virtual ~DeltaAgent();
        /// Cloning an agent.
        virtual Agent* clone() const;
        /// Copy the agent into \c this
        virtual void copy( const Agent & );

        /// Perform one timestep. It consists of checking whether the
        /// agent dies or survives another timestep.
        virtual void step( Population & );
        /// Reproduce. The sibling has its age reset to zero.
        virtual Agent* sibling();
        /// Evaluate network and calculate new attractor
        virtual void evaluate( const Environment & );

        /// Update the state referred to in the evaluation
        void updateReferenceState();
        /// Set reference state
        void referenceState( const boost::dynamic_bitset<> & );
        /// Get reference state
        void referenceState( std::vector< int > & ) const;
        
        /// Write a text representation of the agent to an output stream.
        virtual void writeXml( std::ostream & ) const;

        public:
        /// Set the number of genes for the reference state
        static void nrGenes( int );
        /// Get the number of genes
        static int nrGenes();
        
        protected:
        /// Reference to the final state of the last environment
        boost::dynamic_bitset<> ref_state_;
        
        protected:
        /// Number of genes = max distance of two states in network
        static int nr_genes_;
    };
}
#endif

