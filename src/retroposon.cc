//
// Implementation of a special gene, the Retroposon.
//
// by Anton Crombach, A.B.M.Crombach at uu.nl
//

#include "retroposon.hh"

template<> blowhole::ObjectCache< blowhole::Retroposon >* 
blowhole::ObjectCache< blowhole::Retroposon >::instance_ = 0;


blowhole::ChromosomeElement* 
blowhole::Retroposon::clone() const {
    //return new Retroposon( *this );
    Retroposon *rp = ObjectCache< Retroposon >::instance()->borrowObject();
    rp->copy( *this );
    return rp;
}

void 
blowhole::Retroposon::copy( const ChromosomeElement &ce ) {
    const Retroposon *tp = dynamic_cast< const Retroposon * >( &ce );
    tag_ = tp->tag_;
}

bool
blowhole::Retroposon::toPool() {
    return ObjectCache< Retroposon >::instance()->returnObject( this );
}

void
blowhole::Retroposon::writeXml( std::ostream &os ) const {
    os << "<tposon id=\"" + boost::lexical_cast< std::string >( tag_ ) + "\"/>";
}
