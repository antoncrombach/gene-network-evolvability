//
// Hopfield directed graph
//
// by Anton Crombach, A.B.M.Crombach at uu.nl
//

#ifndef _BLOWHOLE_HOPFIELD_GRAPH_H_
#define _BLOWHOLE_HOPFIELD_GRAPH_H_

#include "defs.hh"
#include "ref_node.hh"

namespace blowhole {

    /// \class HopfieldGraph
    /// \brief Gene regulation network
    class HopfieldGraph : public XmlWriter {
        public:
        /// Typedef edges
        typedef std::pair< int, int > hop_edge;
        /// Typedef of the graph storage structure
        typedef boost::adjacency_list< 
            boost::vecS, boost::vecS, boost::bidirectionalS, 
            boost::property< boost::vertex_name_t, Gene* >,
            boost::property< boost::edge_weight_t, int > > hop_graph;
        typedef boost::adjacency_list< 
            boost::vecS, boost::vecS, boost::bidirectionalS, RefNode, 
            boost::property< boost::edge_weight_t, int > > ref_graph;
        /// Property maps
        /// Hopfield graph
        typedef boost::property_map< hop_graph, boost::vertex_name_t >::type
            name_map;
        typedef boost::property_map< hop_graph, boost::edge_weight_t >::type
            interaction_map;
        /// Accessing vertices and edges
        typedef boost::graph_traits< hop_graph >::vertex_iterator node_iter;
        typedef boost::graph_traits< hop_graph >::in_edge_iterator edge_iter;
        /// and for the reference graph
        typedef boost::graph_traits< ref_graph >::vertex_iterator ref_node_iter;
        typedef boost::graph_traits< ref_graph >::edge_iterator ref_edge_iter;
        // and for the edge mapping
        typedef std::map< hop_edge, int > edge_iaction_map;
        typedef edge_iaction_map::iterator edge_iaction_iter;
        typedef edge_iaction_map::const_iterator edge_iaction_citer;
        
        public:
        /// (Dummy) ctor
        HopfieldGraph();
        /// Copy ctor
        explicit HopfieldGraph( const HopfieldGraph & );
        /// Cloning
        HopfieldGraph* clone() const;
        /// Copying method
        void copy( const HopfieldGraph & );
        /// Dtor
        ~HopfieldGraph();

        /// Build the graph from the genome
        void build( const Genome & );
        /// Calc attractor
        void propagate();
        /// Perturb network
        void perturb();
        /// Hamming distance to other state
        int hamming( const boost::dynamic_bitset<> & ) const;
        /// State as a bitset
        void state( boost::dynamic_bitset<> & ) const;
        
        /// Do we have all nodes?
        bool hasMissingNode() const;
        /// Are we in an attractor?
        bool inAttractor() const;
        /// Are all nodes zero?
        bool allZero() const;
        /// Is network already build?
        bool isEmpty() const;
        
        /// Get graph
        const hop_graph & graph() const;
        
        /// Write in dot format (although it is called xml)
        virtual void writeXml( std::ostream & ) const;

        public:
        /// Noise rate
        static void perturbRate( double );
        /// Set time limit
        static void maxPropagate( int );
        /// Sequential or parallel update
        static void seqPropagate( bool );
        /// Write reference graph to file
        static void writeReferenceGraph( std::ostream & );
        /// Read reference graph from file
        static void readReferenceGraph( std::istream & );
        /// Get reference node ids
        static const std::set< int > & referenceIds();
        /// Get reference nodes
        static const std::map< int, RefNode > & referenceNodes();
        /// Get reference edges
        static const std::map< hop_edge, int > & referenceEdges();
        
        private:
        // Set dynamic_bitset \c state_ to current network state
        void doState();
        // Propagate a single node
        bool propagateNode( hop_graph::vertex_descriptor );
        // Sequential propagating
        void doSeqPropagate();
        // Parallel propagating
        void doParPropagate();

        private:
        // Auxilary function for reading the reference network from file
        static void doReadReferenceGraph( std::istream &is, ref_graph &ref,
            boost::dynamic_properties &dp );
        
        protected:
        mutable hop_graph graph_;
        boost::dynamic_bitset<> state_;
        name_map graph_gmap_;
        interaction_map graph_imap_;
        bool missing_nodes_, in_attractor_;
        
        protected:
        // Time limit for stable state calculation
        static int max_propagate_;
        // Sequential update (true) or parallel (false)
        static bool seq_propagate_;
        // Noise
        static double perturb_rate_;
        // Essential nodes (used for checking if lethal, ie missing a node)
        static std::set< int > ref_ids_;
        // Essential nodes with their info
        static std::map< int, RefNode > ref_nodes_;
        // Which edges have we got (structure of net is predefined)
        static std::map< hop_edge, int > ref_edges_;
    };
    
    inline bool HopfieldGraph::hasMissingNode() const
    { return missing_nodes_; }
    
    inline bool HopfieldGraph::inAttractor() const
    { return in_attractor_; }
    
    inline bool HopfieldGraph::isEmpty() const
    { return boost::num_vertices( graph_ ) == 0; }
    
    inline const HopfieldGraph::hop_graph & HopfieldGraph::graph() const
    { return graph_; }
    
    inline void HopfieldGraph::state( boost::dynamic_bitset<> &v ) const 
    { v = state_; }
}
#endif
