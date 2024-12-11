BEGIN { OFS = "\t" }
{
    prefix = substr($1, 1, 2)
    $4 = prefix $4
    print $0
}
