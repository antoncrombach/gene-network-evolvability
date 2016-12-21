//
// Implementation of abstract agent.
//
// by Anton Crombach, A.B.M.Crombach at uu.nl
//

#include "agent.hh"

blowhole::Agent::Agent() 
: dying_( false ), me_(), ancestor_(), type_( -1 ) {}

blowhole::Agent::Agent( AgentTag t ) 
: dying_( false ), me_( t ), ancestor_(), type_( -1 ) {}

blowhole::Agent::Agent( int tt ) 
: dying_( false ), me_(), ancestor_(), type_( tt ) {}

blowhole::Agent::Agent( const Agent &ag ) {
    copy( ag );
}

void
blowhole::Agent::copy( const Agent &ag ) {
    dying_ = ag.dying_;
    me_ = ag.me_;
    ancestor_ = ag.ancestor_;
    type_ = ag.type_;
}

