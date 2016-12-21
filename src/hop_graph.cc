//
// Implementation of hopfield graph
//
// by Anton Crombach, A.B.M.Crombach at uu.nl
//

#include "hop_graph.hh"
#include "genome.hh"
#include "chromosome.hh"
#include "gene.hh"

int blowhole::HopfieldGraph::max_propagate_ = 0;
bool blowhole::HopfieldGraph::seq_propagate_ = true;
double blowhole::HopfieldGraph::perturb_rate_ = 0.0;
std::set< int > blowhole::HopfieldGraph::ref_ids_;
std::map< int, blowhole::RefNode > blowhole::HopfieldGraph::ref_nodes_;
std::map< blowhole::HopfieldGraph::hop_edge, int >
    blowhole::HopfieldGraph::ref_edges_;

blowhole::HopfieldGraph::HopfieldGraph() 
    : graph_(), state_(), graph_gmap_(), graph_imap_(), 
      missing_nodes_( false ), in_attractor_( false ) {}

blowhole::HopfieldGraph::HopfieldGraph( const HopfieldGraph &hg ) 
    : graph_(), state_(), graph_gmap_(), graph_imap_() {
    copy( hg );
}

blowhole::HopfieldGraph::~HopfieldGraph() {}

blowhole::HopfieldGraph* 
blowhole::HopfieldGraph::clone() const {
    return new HopfieldGraph( *this );
}

void 
blowhole::HopfieldGraph::copy( const HopfieldGraph &hg ) {
    // not copying the graph or the two mappings, as these need to be rebuilt
    state_ = hg.state_;
    missing_nodes_ = hg.missing_nodes_;
    in_attractor_ = hg.in_attractor_;
}

void
blowhole::HopfieldGraph::build( const Genome &g ) {
    // pre: genes and interactions are initialised
    typedef std::multimap< int, hop_graph::vertex_descriptor > vertex_map;
    typedef std::map< Gene*, hop_graph::vertex_descriptor > gene_map;
    
    vertex_map tag_node_map;
    gene_map gene_node_map;
    std::set< int > gene_set;
    // get all genes by id in multimap
    graph_.clear();
    std::list< Chromosome* > ch = g.chromosomes();
    for( Genome::chr_iter ii = ch.begin(); ii != ch.end(); ++ii ) {
        // usually only one chromosome
        std::list< ChromosomeElement* > ces = ( **ii ).elements();
        for( Chromosome::ce_iter jj = ces.begin(); jj != ces.end(); ++jj ) {
            Gene *gg = dynamic_cast< Gene* >( *jj );
            if( gg ) {
                hop_graph::vertex_descriptor aux = 
                    boost::add_vertex( gg, graph_ );
                tag_node_map.insert( std::make_pair( gg->tag(), aux ) );
                gene_node_map.insert( std::make_pair( gg, aux ) );
                gene_set.insert( gg->tag() );
            }
        }
    }

//    cout << g << endl;
//    std::copy( gene_set.begin(), gene_set.end(), 
//        std::ostream_iterator< int >( cout, "\t" ) );
//    cout << endl;
    
    // check if we have all genes
    if( ref_ids_ == gene_set ) {
//        cout << "missing false" << endl;
        missing_nodes_ = false;
        // each gene is preceded by its bindingsites, however now we also have
        // to let repeat elements break up the promotor region...
        // I should not use ref_edges_ here..
        for( Genome::chr_iter ii = ch.begin(); ii != ch.end(); ++ii ) {
            // usually only one chromosome
            Gene* target = 0;
            std::list< ChromosomeElement* > ces = ( **ii ).elements();
            for( Chromosome::ce_riter j = ces.rbegin(); j != ces.rend(); ++j ) {
                if( Chromosome::IsGene()( *j ) ) {
                    target = dynamic_cast< Gene* >( *j );
                } else if( Chromosome::IsInteraction()( *j ) ) {
                    Interaction *ia = dynamic_cast< Interaction* >( *j );
                    vertex_map::iterator a, b, k;
                    boost::tie( a, b ) = tag_node_map.equal_range( ia->tag() );
                    for( k = a; k != b; ++k ) {
                        boost::add_edge( k->second, gene_node_map[ target ], 
                            ia->weight(), graph_ );
                    }
                }
            }
        }
        // cache the gene and interaction maps, not gonna change until rebuilt
        graph_gmap_ = boost::get( boost::vertex_name, graph_ );
        graph_imap_ = boost::get( boost::edge_weight, graph_ );
    } else {
//        cout << "missing true" << endl;
        missing_nodes_ = true;
    }
}

void
blowhole::HopfieldGraph::propagate() {
    // commented out the hardly ever used option of sequential updating    
#ifdef SEQUENTIAL
    if( seq_propagate_ ) { 
        doSeqPropagate();
    } else {
#endif
        doParPropagate();
#ifdef SEQUENTIAL
    }
#endif
    // unconditional update of state
    doState();
}

void
blowhole::HopfieldGraph::doParPropagate() {
    std::map< Gene*, int > config;
    // can even cache these, not?
    // update as long as configurations change (so not keeping track of cycles,
    // but the max_propagate_ shouldn't be a large value anyway)
    int k = 0, l = max_propagate_;
    while( k != l ) {
        // store old configuration in \c config
        node_iter i, end;
        for( boost::tie( i, end ) = boost::vertices( graph_ ); i != end; ++i ) {
            Gene *aux = boost::get( graph_gmap_, *i );
            config[ aux ] = aux->state();
        }
        // update nodes in graph according to state in copy
        // and check if any of them really changes (so we may stop)
        bool any_change = false;
        for( boost::tie( i, end ) = boost::vertices( graph_ ); i != end; ++i ) {
            double sum = 0.0;
            edge_iter j, ant;
            for( boost::tie( j, ant ) = boost::in_edges( *i, graph_ ); 
                j != ant; ++j ) {
                Gene *aux = boost::get( graph_gmap_, 
                    boost::source( *j, graph_ ) );
                sum += config[ aux ] * boost::get( graph_imap_, *j );
            }
            Gene *aux = boost::get( graph_gmap_, *i );
            // boolean as int gives: if true then 1 else 0
            /*if( !close_to( sum, 0.0 ) ) {*/
            if( !close_to( sum, static_cast< double >( aux->threshold() ) ) ) {
                aux->state( static_cast< int >( sum > aux->threshold() ) );
                any_change |= ( config[ aux ] != aux->state() );
            }
        }
        // proceed or stop...
        if( any_change ) {
            in_attractor_ = false;
            ++k;
        } else {
            in_attractor_ = true;
            l = k;
        }
    }
}

void
blowhole::HopfieldGraph::doSeqPropagate() {
#ifdef SEQUENTIAL
    const int nn = boost::num_vertices( graph_ );
    std::vector< hop_graph::vertex_descriptor > aux;
    aux.reserve( nn );
    node_iter i, end;
    for( boost::tie( i, end ) = boost::vertices( graph_ ); i != end; ++i ) {
        aux.push_back( *i );
    }
    std::random_shuffle( aux.begin(), aux.end(), rand_range< int > );
    int j = 0, k = 0;
    int last = 0;
    do {
        if( k == nn ) {
            std::random_shuffle( aux.begin(), aux.end(), rand_range< int > );
            k = 0;
        }
        ++k;
        ++j;        
        if( propagateNode( aux[ k ] ) ) {
            last = j;
        }
    } while( j - last < static_cast< int >( 10 * nn ) and j < max_propagate_ );
    // flag if we didn't hit a stable state?
    in_attractor_ = false;
#endif
}

bool
blowhole::HopfieldGraph::propagateNode( hop_graph::vertex_descriptor v ) {
#ifdef SEQUENTIAL
    // initial value of sum matters...
    double sum = 0.0;
    edge_iter i, end;
    for( boost::tie( i, end ) = boost::in_edges( v, graph_ ); i != end; ++i ) {
        Gene *aux = boost::get( graph_gmap_, boost::source( *i, graph_ ) );
        sum += aux->state() * boost::get( graph_imap_, *i );
    }
    Gene *aux = boost::get( graph_gmap_, v );
    // boolean as int gives: if true then 1 else 0
    int bux = aux->state();
    // update according to Li et. al. "The yeast cell-cycle network
    // is robustly designed"
    /*if( !close_to( sum, 0.0 ) ) {*/
    if( !close_to( sum, static_cast< double >( aux->threshold() ) ) ) {
        aux->state( static_cast< int >( sum > aux->threshold() ) );
    }
    return bux != aux->state();
#else
    return false;
#endif
}

int
blowhole::HopfieldGraph::hamming( const boost::dynamic_bitset<> &b ) const {
    // binary hamming distance between two network states (not configurations!)
    boost::dynamic_bitset<> aux( state_ );    
    // note: sizes should be equal.. if not, count the ones in the extra part
    // of our state (attractor state is zero in this region)
    int extra = 0;
    if( b.size() < aux.size() ) {
        for( uint i = b.size(); i != aux.size(); ++i ) {
            if( aux.test( i ) ) ++extra;
        }
    }
    aux.resize( b.size() );
    aux ^= b;
    return aux.count() + extra;
}

void
blowhole::HopfieldGraph::perturb() {
    if( close_to( perturb_rate_, 0.0 ) ) return;
    node_iter i, end;
    for( boost::tie( i, end ) = boost::vertices( graph_ ); i != end; ++i ) {
        Gene *aux = boost::get( graph_gmap_, *i );
        if( uniform() < perturb_rate_ ) {
            aux->mutateState();
            in_attractor_ = false;
        }
    }
}

void
blowhole::HopfieldGraph::doState() {
    // pre: vertex ids element of [0..n)
    std::vector< double > bux( ref_ids_.size(), 0.0 );
    node_iter i, end;
    for( boost::tie( i, end ) = boost::vertices( graph_ ); i != end; ++i ) {
        Gene *aux = boost::get( graph_gmap_, *i );
        bux[ aux->tag() ] += aux->state();
    }
    state_.resize( ref_ids_.size() );
    for( uint jj = 0; jj < bux.size(); ++jj ) {
        // virtual genes behave differently...
        // "at least one gene is on"
        state_[ jj ] = ( bux[ jj ] > 0.5 );
    }
}

bool
blowhole::HopfieldGraph::allZero() const {
    return state_.none();
}    

void
blowhole::HopfieldGraph::writeXml( std::ostream &os ) const {
    // NOTE: be aware that the id's of this graph do not match the reference
    // graph, it is a permutation...
    os << "<network>\n";
    boost::dynamic_properties dp;
    dp.property( "node_id", boost::get( boost::vertex_index, graph_ ) );
    dp.property( "iaction", boost::get( boost::edge_weight, graph_ ) );
    boost::write_graphviz( os, graph_, dp );
    os << "</network>\n";
}

void
blowhole::HopfieldGraph::perturbRate( double t ) {
    perturb_rate_ = t;
}

void
blowhole::HopfieldGraph::maxPropagate( int t ) {
    max_propagate_ = t;
}

void
blowhole::HopfieldGraph::seqPropagate( bool t ) {
    seq_propagate_ = t;
}

void
blowhole::HopfieldGraph::readReferenceGraph( std::istream &is ) {
    ref_graph ref;
    boost::dynamic_properties dp;
    doReadReferenceGraph( is, ref, dp );
    
    // now create set of reference ids and nodes
    {
        ref_node_iter ii, end;
        for( boost::tie( ii, end ) = boost::vertices( ref ); ii != end; ++ii ) {
            // get node id and store it
            int aux = boost::get( &RefNode::node_id, ref, *ii );
            ref_ids_.insert( aux );
            ref_nodes_.insert( std::make_pair( 
                aux, boost::get( boost::vertex_bundle, ref, *ii ) ) );
        }
    }
    // and set of reference edges
    {
        ref_edge_iter ii, end;
        for( boost::tie( ii, end ) = boost::edges( ref ); ii != end; ++ii ) {
            // ref_edge is mapping from edge to interaction
            ref_edges_[ std::make_pair(
                boost::source( *ii, ref ), boost::target( *ii, ref ) ) ] = 
                boost::get( boost::edge_weight, ref, *ii );
        }
    }
}

const std::set< int > &
blowhole::HopfieldGraph::referenceIds() {
    return ref_ids_;
}

const std::map< int, blowhole::RefNode > &
blowhole::HopfieldGraph::referenceNodes() {
    return ref_nodes_;
}

const std::map< blowhole::HopfieldGraph::hop_edge, int > &
blowhole::HopfieldGraph::referenceEdges() {
    return ref_edges_;
}

void
blowhole::HopfieldGraph::writeReferenceGraph( std::ostream &os ) {
    typedef std::map< int, ref_graph::vertex_descriptor > vertex_map;
    
    ref_graph ref;
    vertex_map gene_node_map;
    // build graph from ref_nodes_ and ref_edges_
    for( std::map< int, RefNode >::const_iterator ii = ref_nodes_.begin();
        ii != ref_nodes_.end(); ++ii ) {
        gene_node_map.insert( 
            std::make_pair( ii->first, boost::add_vertex( ii->second, ref ) ) );
    }
    // and the edges
    for( edge_iaction_citer i = ref_edges_.begin(); i != ref_edges_.end(); ++i){
        boost::add_edge( gene_node_map[ i->first.first ],
            gene_node_map[ i->first.second ], i->second, ref );
    }

    boost::dynamic_properties dp;
    dp.property( "node_id", boost::get( &RefNode::node_id, ref ) );
    dp.property( "thr", boost::get( &RefNode::threshold, ref ) );
    dp.property( "iaction", boost::get( boost::edge_weight, ref ) );
    boost::write_graphviz( os, ref, dp );
}    

void
blowhole::HopfieldGraph::doReadReferenceGraph( std::istream &is, ref_graph &ref,
    boost::dynamic_properties &dp ) {
    // auxilary function
    dp.property( "node_id", boost::get( &RefNode::node_id, ref ) );
    dp.property( "thr", boost::get( &RefNode::threshold, ref ) );
    dp.property( "iaction", boost::get( boost::edge_weight, ref ) );

    try {
        boost::read_graphviz( is, ref, dp );
    } catch( const boost::graph_exception &e ) {
        std::cout << "Exception: " << e.what() << "\n";
    } catch( const boost::dynamic_property_exception &e ) {
        std::cout << "Exception: " << e.what() << "\n";
    }
}


// Leftovers
/*double
blowhole::HopfieldGraph::hamming( const std::vector< double > &v ) const {
    std::vector< double > &s;
    stateAsVector( s );
    // if size different...
    std::vector< double >::iterator min;
    int aux = s.size() - v.size();
    double extra = 0.0;
    if( aux > 0 ) {
        // padding zeroes to s (well, just count the extra difference)
        min = boost::next( s.begin(), v.size() );
        extra = std::accumulate( min, s.end(), 0.0 );
    } else if( aux < 0 ) {
        // padding zeroes to v (same trick)
        min = boost::next( v.begin(), s.size() );
        extra = std::accumulate( min, v.end(), 0.0 );
        min = s.end();
    }
    // subtract, absolute values and accumulate difference
    // reuse iterator \c i as end
    std::transform( s.begin(), min, v.begin(), s.begin(), 
        std::compose1( std::minus< double >(), absolute< double >() ) );
    return std::accumulate( s.begin(), min, extra );
}*/

/*void
blowhole::HopfieldGraph::state( std::vector< double > &v ) {
    // pre: v.empty
    for( boost::dynamic_bitset<>::size_type i = 0; i < state_.size(); ++i ) {
        v.push_back( state_[ i ] );
    }
}*/

