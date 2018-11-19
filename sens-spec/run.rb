#!/usr/bin/env ruby
require 'optparse'
require 'fileutils'
require '/home/cew54/slurmer/Qsub.rb'
require 'pp'

## empty array for commands

OPTIONS = {}
OPTIONS[:int] = false
OPTIONS[:autorun] = false
OPTIONS[:all] = false
OptionParser.new do |opts|
  opts.banner = "Usage: runguess.rb [OPTIONS] COMMAND [DIR]"

  # opts.on("-a", "--[no-]all", "(Re-)run all outputs regardless of whether it already exists") do |i|
  #   OPTIONS[:all] = i
  # end
  
  opts.on("-i", "--[no-]interactive", "Run each job in serial on the interactive log in node") do |i|
    OPTIONS[:int] = i
  end
  
  opts.on("-r", "--autoRun", "Run each job on the queue without user input") do |i|
    OPTIONS[:autorun] = i
  end
  
  opts.on("-v", "--[no-]verbose", "Run verbosely") do |v|
    OPTIONS[:verbose] = v
  end

  opts.on("-h", "--help", "Show help") do |h|
    OPTIONS[:help] = h
  end
end.parse!
COMMAND = ARGV.shift


 def usage()
  puts OPTIONS
  puts "Usage: run.rb [OPTIONS] COMMAND

  COMMANDS are:
      run   :: run additional reps as needed
      summ  :: summarise current state
      clean :: remove slurm files      
  "
end

 if OPTIONS[:help] then
  usage()
  exit 0
end
commands = []
COMFILE="tmprun.sh"
args = {:job=>COMMAND,
        :time=>"12:00:00",
        :tasks=>'1',
        :cpus=>'2',
        :autorun=>OPTIONS[:autorun],
        :excl=>" "}

## all baits todo
DIR="/home/cew54/scratch/peaky"
allbaits = File.readlines("#{DIR}/totest.txt").map(&:chomp)

PATTERNS=%w{ N0 N1 N2 N0.N1 N1.N2 N0.N2 }

PATTERNS.each {|p|

  done = Dir.glob("#{DIR}/#{p}/*.csv").map{ |x| File.basename(x,".csv") }
  todo = allbaits - done
  n=todo.length()
  of="#{DIR}/#{p}-totest.txt"
  File.open(of, "w") do |f|
    f.puts(todo)
  end
  puts "pattern #{p}, todo = #{n}"
  flag = OPTIONS[:autorun] ? "-r" : ""
  comm = "qR.rb #{flag} -c 1 -t \"8:00:00\" -j #{p} -p skylake-himem -a CWALLACE-SL2-CPU -y 1-#{n} ./do-tests.R --args patt=#{p} file=#{of}"
  system( comm )

}
