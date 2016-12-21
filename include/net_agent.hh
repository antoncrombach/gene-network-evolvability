//
// Interface of an simple agent.
//
// by Anton Crombach, A.B.M.Crombach at uu.nl
//

#ifndef _BLOWHOLE_NETAGENT_H_
#define _BLOWHOLE_NETAGENT_H_

#include "defs.hh"
#include "genreg_agent.hh"
#include "genome.hh"
#include "hop_graph.hh"

namespace blowhole {

    /// \class NetAgent
    /// \brief Simple agent with a network.
    class NetAgent : public GenRegAgent {
        public:
        /// Default ctor.
        NetAgent();
        /// Ctor.
        NetAgent( uint tag, Genome *g );
        /// Copy ctor.
        NetAgent( const NetAgent & );
        /// Destructor.
        virtual ~NetAgent();
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
        
        /// Write a text representation of the agent to an output stream.
        virtual void writeXml( std::ostream & ) const;
    };
}
#endif

